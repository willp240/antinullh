// Antinu headers
#include <PDFConfigLoader.hh>
#include <DistBuilder.hh>
#include <FitConfigLoader.hh>
#include <EventConfigLoader.hh>
#include <SystConfigLoader.hh>
#include <SystFactory.hh>
#include <OscGridConfigLoader.hh>

// OXO headers
#include <ROOTNtuple.h>
#include <BinnedNLLH.h>
#include <IO.h>
#include <Rand.h>

// ROOT headers
#include <TH1D.h>

// c++ headers
#include <sys/stat.h>

using namespace antinufit;

void llh_scan(const std::string &fitConfigFile_,
              const std::string &evConfigFile_,
              const std::string &pdfConfigFile_,
              const std::string &systConfigFile_,
              const std::string &oscGridConfigFile_)
{
  Rand::SetSeed(0);

  // Load up the fit configuration information
  FitConfig fitConfig;
  FitConfigLoader fitLoader(fitConfigFile_);
  fitConfig = fitLoader.LoadActive();
  bool isAsimov = fitConfig.GetAsimov();
  double livetime = fitConfig.GetLivetime();
  bool beestonBarlowFlag = fitConfig.GetBeestonBarlow();
  std::string outDir = fitConfig.GetOutDir();
  ParameterDict constrMeans = fitConfig.GetConstrMeans();
  ParameterDict constrSigmas = fitConfig.GetConstrSigmas();
  ParameterDict mins = fitConfig.GetMinima();
  ParameterDict maxs = fitConfig.GetMaxima();
  ParameterDict noms = fitConfig.GetNominals();

  struct stat st = {0};
  if (stat(outDir.c_str(), &st) == -1)
    mkdir(outDir.c_str(), 0700);

  // Load up all the event types we want to contribute
  typedef std::map<std::string, EventConfig> EvMap;
  EventConfigLoader loader(evConfigFile_);
  EvMap toGet = loader.LoadActive();

  // Load up the PDF information (skeleton axis details, rather than the distributions themselves)
  PDFConfigLoader pdfLoader(pdfConfigFile_);
  PDFConfig pdfConfig = pdfLoader.Load();
  std::string pdfDir = pdfConfig.GetPDFDir();
  std::vector<std::string> dataObs = pdfConfig.GetDataBranchNames();
  ObsSet dataObsSet(dataObs);
  AxisCollection systAxes = DistBuilder::BuildAxes(pdfConfig, dataObs.size());

  // Load up the systematics
  SystConfigLoader systLoader(systConfigFile_);
  SystConfig systConfig = systLoader.LoadActive();
  std::map<std::string, std::string> systParamNames = systConfig.GetParamNames();
  std::map<std::string, std::string> systGroup = systConfig.GetGroup();
  std::map<std::string, std::string> systType = systConfig.GetType();
  std::map<std::string, std::vector<std::string>> systDistObs = systConfig.GetDistObs();
  std::map<std::string, std::vector<std::string>> systTransObs = systConfig.GetTransObs();
  std::vector<std::string> fullParamNameVec;

  // Load up the oscillation probability grids
  OscGridConfigLoader oscGridLoader(oscGridConfigFile_);
  OscGridConfig oscGridConfig = oscGridLoader.Load();
  std::string outfilename = oscGridConfig.GetFilename();
  std::string reactorjson = oscGridConfig.GetReactorsJsonFile();
  std::map<int, OscGrid *> oscGridMap;

  // First read the reactor distance info
  std::unordered_map<int, double> indexDistance = LoadIndexDistanceMap(reactorjson);
  for (std::unordered_map<int, double>::iterator it = indexDistance.begin(); it != indexDistance.end(); ++it)
  {
    std::string oscGridFileName = outfilename + "_" + std::to_string(it->first) + ".root";
    OscGrid *oscGrid = new OscGrid(oscGridFileName);
    oscGridMap[it->first] = oscGrid;
  }

  // Loop over systematics and declare each one
  std::map<std::string, Systematic *> systMap;
  for (std::map<std::string, std::string>::iterator it = systType.begin(); it != systType.end(); ++it)
  {
    std::cout << "Constructing " << it->first << std::endl;
    std::vector<std::string> paramNameVec;
    std::stringstream ss(systParamNames[it->first]);
    std::string paramName;
    while (std::getline(ss, paramName, ','))
    {
      paramNameVec.push_back(paramName);
      if (noms.find(paramName) == noms.end())
      {
        std::cout << "ERROR: Syst config defines systematic with parameter " << paramName << ", but this parameter is not defined in the fit config" << std::endl;
        throw;
      }
    }
    fullParamNameVec.insert(fullParamNameVec.end(), paramNameVec.begin(), paramNameVec.end());
    Systematic *syst = SystFactory::New(it->first, systType[it->first], paramNameVec, noms, oscGridMap, indexDistance);
    AxisCollection systAxes = DistBuilder::BuildAxes(pdfConfig, systDistObs[it->first].size());
    syst->SetAxes(systAxes);
    // The "dimensions" the systematic applies to
    syst->SetTransformationObs(systTransObs[it->first]);
    // All the "dimensions" of the dataset
    syst->SetDistributionObs(systDistObs[it->first]);
    syst->Construct();
    systMap[it->first] = syst;
  }

  // A parameter could have been defined in the fit config but isn't associated with a pdf or systematic
  // If so ignore that parameter
  ParameterDict::iterator it = noms.begin();
  size_t numParams = noms.size();
  for (int iParam = 0; iParam < numParams; iParam++)
  {
    bool iterate = true;
    if (toGet.find(it->first) == toGet.end() && std::find(fullParamNameVec.begin(), fullParamNameVec.end(), it->first) == fullParamNameVec.end())
    {
      std::cout << it->first << " parameter defined in fit config but not in syst or event config. It will be ignored." << std::endl;
      constrSigmas.erase(it->first);
      constrMeans.erase(it->first);
      mins.erase(it->first);
      maxs.erase(it->first);
      it = noms.erase(it);
      if (it != noms.begin())
        it--;
      else
        iterate = false;
    }
    if (iterate)
      it++;
  }

  // Create the individual PDFs and Asimov components (could these be maps not vectors?)
  std::vector<BinnedED> pdfs;
  std::vector<BinnedED> indivAsmvDists;
  std::vector<int> genRates;
  std::vector<std::vector<std::string>> pdfGroups;

  // Create the empty full dist
  BinnedED asimov = BinnedED("asimov", systAxes);
  asimov.SetObservables(dataObs);

  // Build each of the PDFs, scale them to the correct size
  for (EvMap::iterator it = toGet.begin(); it != toGet.end(); ++it)
  {

    std::cout << "Building distribution for " << it->first << std::endl;
    pdfGroups.push_back(it->second.GetGroup());

    // Get MC events from file
    DataSet *dataSet;
    try
    {
      std::cout << "Loading " <<  it->second.GetPrunedPath() << std::endl;
      dataSet = new ROOTNtuple(it->second.GetPrunedPath(), "pruned");
    }
    catch (const IOError &e_)
    {
      std::cout << "Warning: skipping " << it->first << " couldn't open files:\n\t" << e_.what() << std::endl;
      continue;
    }

    // Build distribution of those events
    BinnedED dist;
    int num_dimensions = it->second.GetNumDimensions();
    dist = DistBuilder::Build(it->first, num_dimensions, pdfConfig, dataSet);
    dist.AddPadding(1E-6);

    // Save the generated number of events for Beeston Barlow
    genRates.push_back(dist.Integral());

    // Scale for PDF and add to vector
    if (dist.Integral())
      dist.Normalise();
    pdfs.push_back(dist);

    // Apply nominal systematic variables
    for (std::map<std::string, Systematic *>::iterator it = systMap.begin(); it != systMap.end(); ++it)
    {
      // If group is "", we apply to all groups
      if (systGroup[it->first] == "" || std::find(pdfGroups.back().begin(), pdfGroups.back().end(), systGroup[it->first]) != pdfGroups.back().end())
      {
        double distInt = dist.Integral();
        dist = it->second->operator()(dist);
        dist.Scale(distInt);
      }
    }

    // Now scale the Asimov component by expected count, and also save pdf as a histo
    dist.Scale(noms[it->first]);
    if (dist.GetNDims() != asimov.GetNDims())
    {
      BinnedED marginalised = dist.Marginalise(dataObs);
      asimov.Add(marginalised);
      indivAsmvDists.push_back(marginalised);
      IO::SaveHistogram(marginalised.GetHistogram(), pdfDir + "/" + it->first + ".root", dist.GetName());
    }
    else
    {
      asimov.Add(dist);
      indivAsmvDists.push_back(dist);
      IO::SaveHistogram(dist.GetHistogram(), pdfDir + "/" + it->first + ".root", dist.GetName());
    }
  }

  // And save combined histogram (final asimov dataset)
  IO::SaveHistogram(asimov.GetHistogram(), outDir + "/asimov.root", "asimov");

  // Now let's load up the data
  BinnedED dataDist;
  if (!isAsimov)
  {

    std::string dataPath = fitConfig.GetDatafile();

    // Could be h5 or root file
    if (dataPath.substr(dataPath.find_last_of(".") + 1) == "h5")
    {
      Histogram loaded = IO::LoadHistogram(dataPath);
      dataDist = BinnedED("data", loaded);
      dataDist.SetObservables(pdfConfig.GetBranchNames());
    }
    else
    {
      // Load up the data set
      ROOTNtuple dataToFit(dataPath, "pruned");

      // And bin the data inside
      dataDist = DistBuilder::Build("data", pdfConfig, (DataSet *)&dataToFit);
    }
  }
  else
    dataDist = asimov;

  // Now build the likelihood
  BinnedNLLH lh;
  // Set the buffer region
  // lh.SetBufferAsOverflow(true);
  //lh.SetBuffer("energy", 10, 10);
  // Add our PDFs and data
  lh.AddPdfs(pdfs, pdfGroups, genRates);
  lh.SetDataDist(dataDist);
  // Set whether or not to use Beeston Barlow
  lh.SetBarlowBeeston(beestonBarlowFlag);
  // Add the systematics and any prior constraints
  for (std::map<std::string, Systematic *>::iterator it = systMap.begin(); it != systMap.end(); ++it)
  {
    if (systGroup[it->first] != "")
      lh.AddSystematic(it->second, systGroup[it->first]);
    else
      lh.AddSystematic(it->second);
  }

  for (ParameterDict::iterator it = constrMeans.begin(); it != constrMeans.end(); ++it)
    lh.SetConstraint(it->first, it->second, constrSigmas.at(it->first));

  // And finally bring it all together
  lh.RegisterFitComponents();

  // Now onto the LLH Scan. First set the number of points in the scan
  int npoints = 150;
  int countwidth = double(npoints) / double(5);

  // Initialise to nominal values
  ParameterDict parameterValues;
  for (ParameterDict::iterator it = mins.begin(); it != mins.end(); ++it)
    parameterValues[it->first] = noms[it->first];

  // If we have constraints, initialise to those
  for (ParameterDict::iterator it = constrMeans.begin(); it != constrMeans.end(); ++it)
    parameterValues[it->first] = constrMeans[it->first];

  // Set to these initial values
  lh.SetParameters(parameterValues);

  // Setup outfile
  std::string outFileName = outDir + "/llh_scan.root";
  TFile *outFile = new TFile(outFileName.c_str(), "recreate");

  // Loop over parameters
  for (ParameterDict::iterator it = mins.begin(); it != mins.end(); ++it)
  {
    // Get param name
    std::string name = it->first;
    std::cout << "Scanning for " << name << std::endl;

    // Set scan range from this parameter's max and min values
    // We plot x axis as relative to the nominal value, so if this = 0 we set it to 1
    double nom = noms[name];
    if (nom == 0)
      nom = 1;
    double min = mins[name];
    double max = maxs[name];
    // Make histogram for this parameter
    TString htitle = Form("%s, Asimov Rate: %f", name.c_str(), nom);
    TH1D *hScan = new TH1D((name + "_full").c_str(), (name + "_full").c_str(), npoints, min / nom, max / nom);
    hScan->SetTitle(std::string(htitle + ";" + name + " (rel. to Asimov); -(ln L_{full})").c_str());

    // Commented out bits here are in anticipation of splitting LLH calculation in oxo at some point. Could then see prior and sample contributions separately
    // TH1D *hScanSam = new TH1D((name+"_sam").c_str(), (name+"_sam").c_str(), npoints, min/nom, max/nom);
    // hScanSam->SetTitle(std::string(std::string("2LLH_sam, ") + name + ";" + name + "; -2(ln L_{sample})").c_str());
    // TH1D *hScanPen = new TH1D((name+"_pen").c_str(), (name+"_pen").c_str(), npoints, min/nom, max/nom);
    // hScanPen->SetTitle(std::string(std::string("2LLH_pen, ") + name + ";" + name + "; -2(ln L_{penalty})").c_str());

    // Now loop from min to max in npoint steps
    for (int i = 0; i < npoints; i++)
    {

      if (i % countwidth == 0)
        std::cout << i << "/" << npoints << " (" << double(i) / double(npoints) * 100 << "%)" << std::endl;

      // Set Parameters
      double parval = hScan->GetBinCenter(i + 1) * nom;
      double tempval = parameterValues[name];
      parameterValues[name] = parval;
      lh.SetParameters(parameterValues);

      // Eval LLH (later do sample and penalty)
      double llh = lh.Evaluate();

      // SetBinContents
      hScan->SetBinContent(i + 1, llh);
      // hScanSam->SetBinContent(i+1, llh);
      // hScanPen->SetBinContent(i+1, llh);

      // return to nominal value
      parameterValues[name] = tempval;
    }
    // Write Histos
    hScan->Write();
    // hScanSam->Write();
    // hScanPen->Write();
  }
  // Close file
  outFile->Close();
  std::cout << "Wrote scan to " << outFile->GetName() << std::endl;

  delete outFile;
}

int main(int argc, char *argv[])
{
  if (argc != 6)
  {
    std::cout << "\nUsage: llh_scan <fit_config_file> <eve_config_file> <pdf_config_file> <syst_config_file> <oscgrid_config_file>" << std::endl;
    return 1;
  }

  std::string fitConfigFile(argv[1]);
  std::string eveConfigFile(argv[2]);
  std::string pdfConfigPath(argv[3]);
  std::string systConfigFile(argv[4]);
  std::string oscgridConfigFile(argv[5]);

  llh_scan(fitConfigFile, eveConfigFile, pdfConfigPath, systConfigFile, oscgridConfigFile);

  return 0;
}
