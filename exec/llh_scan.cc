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
  bool isFakeData = fitConfig.GetFakeData();
  bool beestonBarlowFlag = fitConfig.GetBeestonBarlow();
  std::string outDir = fitConfig.GetOutDir();
  ParameterDict constrMeans = fitConfig.GetConstrMeans();
  ParameterDict constrSigmas = fitConfig.GetConstrSigmas();
  ParameterDict mins = fitConfig.GetMinima();
  ParameterDict maxs = fitConfig.GetMaxima();
  ParameterDict noms = fitConfig.GetNominals();
  ParameterDict constrRatioMeans = fitConfig.GetConstrRatioMeans();
  ParameterDict constrRatioSigmas = fitConfig.GetConstrRatioSigmas();
  std::map<std::string, std::string> constrRatioParName = fitConfig.GetConstrRatioParName();
  ParameterDict constrCorrs = fitConfig.GetConstrCorrs();
  std::map<std::string, std::string> constrCorrParName = fitConfig.GetConstrCorrParName();
  ParameterDict fdValues = fitConfig.GetFakeDataVals();

  std::string pdfDir = outDir + "/unscaled_pdfs";
  std::string asimovDistDir = outDir + "/asimov_dists";
  std::string fakedataDistDir = outDir + "/fakedata_dists";
  // Create output directories
  struct stat st = {0};
  if (stat(outDir.c_str(), &st) == -1)
    mkdir(outDir.c_str(), 0700);
  if (stat(pdfDir.c_str(), &st) == -1)
    mkdir(pdfDir.c_str(), 0700);
  if (stat(asimovDistDir.c_str(), &st) == -1)
    mkdir(asimovDistDir.c_str(), 0700);
  if (stat(fakedataDistDir.c_str(), &st) == -1)
    mkdir(fakedataDistDir.c_str(), 0700);

  // Load up all the event types we want to contribute
  typedef std::map<std::string, EventConfig> EvMap;
  EventConfigLoader loader(evConfigFile_);
  EvMap toGet = loader.LoadActive();

  // Load up the PDF information (skeleton axis details, rather than the distributions themselves)
  PDFConfigLoader pdfLoader(pdfConfigFile_);
  PDFConfig pdfConfig = pdfLoader.Load();
  std::vector<std::string> dataObs = pdfConfig.GetDataBranchNames();
  ObsSet dataObsSet(dataObs);
  AxisCollection systAxes = DistBuilder::BuildAxes(pdfConfig, dataObs.size());

  // Load up the systematics
  SystConfigLoader systLoader(systConfigFile_);
  SystConfig systConfig = systLoader.LoadActive();
  std::map<std::string, std::vector<std::string>> systParamNames = systConfig.GetParamNames();
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
    // Check any parameter defined for a systematic, was declared in the fit config
    for (std::map<std::string, std::vector<std::string>>::iterator paramIt = systParamNames.begin(); paramIt != systParamNames.end(); ++paramIt)
    {
      for (int iParam = 0; iParam < paramIt->second.size(); iParam++)
      {
        if (noms.find(paramIt->second.at(iParam)) == noms.end())
        {
          std::cout << "ERROR: Syst config defines systematic with parameter " << paramIt->second.at(iParam) << ", but this parameter is not defined in the fit config" << std::endl;
          throw;
        }
      }
    }
    fullParamNameVec.insert(fullParamNameVec.end(), systParamNames[it->first].begin(), systParamNames[it->first].end());

    // Now build the systematic
    Systematic *syst = SystFactory::New(it->first, systType[it->first], systParamNames[it->first], noms, oscGridMap, indexDistance);
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

  PrintParams(noms, mins, maxs, constrMeans, constrSigmas, constrRatioMeans, constrRatioSigmas, constrRatioParName, constrCorrs, constrCorrParName);

  // Create the individual PDFs and Asimov components
  std::vector<BinnedED> pdfs;
  std::vector<int> genRates;
  std::vector<std::vector<std::string>> pdfGroups;
  std::vector<NormFittingStatus> *norm_fitting_statuses = new std::vector<NormFittingStatus>;

  // Create the empty full dist
  BinnedED asimov = BinnedED("asimov", systAxes);
  asimov.SetObservables(dataObs);
  // And an empty fake data dist
  BinnedED fakeDataset = BinnedED("fake_dataset", systAxes);
  fakeDataset.SetObservables(dataObs);

  // Build each of the PDFs, scale them to the correct size
  for (EvMap::iterator it = toGet.begin(); it != toGet.end(); ++it)
  {

    std::cout << "Building distribution for " << it->first << std::endl;
    pdfGroups.push_back(it->second.GetGroup());

    // Get MC events from file
    DataSet *dataSet;
    try
    {
      std::cout << "Loading " << it->second.GetPrunedPath() << std::endl;
      dataSet = new ROOTNtuple(it->second.GetPrunedPath(), "pruned");
    }
    catch (const IOError &e_)
    {
      std::cout << "Warning: skipping " << it->first << " couldn't open files:\n\t" << e_.what() << std::endl;
      continue;
    }

    // Build distribution of those events
    BinnedED dist;
    BinnedED fakeDataDist;
    BinnedED unscaledPDF;
    int num_dimensions = it->second.GetNumDimensions();
    dist = DistBuilder::Build(it->first, num_dimensions, pdfConfig, dataSet);

    // Now make a fake data dist for the event type
    fakeDataDist = dist;

    // Save the generated number of events for Beeston Barlow
    genRates.push_back(dist.Integral());

    // Scale for PDF and add to vector
    if (dist.Integral() && fakeDataDist.Integral())
    {
      dist.Normalise();
      fakeDataDist.Normalise();
    }
    pdfs.push_back(dist);
    unscaledPDF = dist;

    norm_fitting_statuses->push_back(INDIRECT);

    // Now make fake data dist for the event type
    fakeDataDist = dist;

    double distInt = dist.Integral();
    double norm;
    // Apply nominal systematic variables
    for (std::map<std::string, Systematic *>::iterator systIt = systMap.begin(); systIt != systMap.end(); ++systIt)
    {
      // If group is "", we apply to all groups
      if (systGroup[systIt->first] == "" || std::find(pdfGroups.back().begin(), pdfGroups.back().end(), systGroup[systIt->first]) != pdfGroups.back().end())
      {
        dist = systIt->second->operator()(dist, &norm);
        // Set syst parameter to fake data value, and apply to fake data dist and rescale
        std::set<std::string> systParamNames = systMap[systIt->first]->GetParameterNames();
        for (auto itSystParam = systParamNames.begin(); itSystParam != systParamNames.end(); ++itSystParam)
        {
          systMap[systIt->first]->SetParameter(*itSystParam, fdValues[*itSystParam]);
        }
        systIt->second->Construct();
        fakeDataDist = systIt->second->operator()(fakeDataDist, &norm);
      }
    }
    dist.Scale(distInt);
    fakeDataDist.Scale(distInt);

    // Now scale the Asimov component by expected count, and also save pdf as a histo
    dist.Scale(noms[it->first]);
    fakeDataDist.Scale(fdValues[it->first]);
    if (dist.GetNDims() != asimov.GetNDims())
    {
      BinnedED marginalised = dist.Marginalise(dataObs);
      asimov.Add(marginalised);
      BinnedED marginalisedPDF = unscaledPDF.Marginalise(dataObs);
      IO::SaveHistogram(marginalised.GetHistogram(), asimovDistDir + "/" + it->first + ".root", dist.GetName());
      IO::SaveHistogram(marginalisedPDF.GetHistogram(), pdfDir + "/" + it->first + ".root", unscaledPDF.GetName());

      // Also scale fake data dist by fake data value
      if (isFakeData)
      {
        BinnedED marginalisedFakeData = fakeDataDist.Marginalise(dataObs);
        fakeDataset.Add(marginalisedFakeData);
        IO::SaveHistogram(marginalisedFakeData.GetHistogram(), fakedataDistDir + "/" + it->first + ".root", fakeDataDist.GetName());
      }
    }
    else
    {
      asimov.Add(dist);
      IO::SaveHistogram(dist.GetHistogram(), asimovDistDir + "/" + it->first + ".root", dist.GetName());
      IO::SaveHistogram(unscaledPDF.GetHistogram(), pdfDir + "/" + it->first + ".root", unscaledPDF.GetName());
      // Also scale fake data dist by fake data value
      if (isFakeData)
      {
        fakeDataset.Add(fakeDataDist);
        IO::SaveHistogram(fakeDataDist.GetHistogram(), fakedataDistDir + "/" + it->first + ".root", fakeDataDist.GetName());
      }
    }
  }

  // And save combined histogram (final asimov dataset). This is with everything (including oscillation parameters) nominal. This is what we
  // compare to to calculate the llh at each point in the scans
  IO::SaveHistogram(asimov.GetHistogram(), outDir + "/asimov.root", "asimov");
  if (isFakeData)
  {
    IO::SaveHistogram(fakeDataset.GetHistogram(), outDir + "/fakedata.root", "fakedata");
  }

  // Now let's load up the data
  BinnedED dataDist;
  if (isAsimov)
  {
    dataDist = asimov;
  }
  else if (isFakeData)
  {
    dataDist = fakeDataset;
  }
  else
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

  // Now build the likelihood
  BinnedNLLH lh;
  lh.SetBuffer("energy", 1, 20);
  // Add our data
  lh.SetDataDist(dataDist);
  // Set whether or not to use Beeston Barlow
  lh.SetBarlowBeeston(beestonBarlowFlag);
  // Add the systematics
  for (std::map<std::string, Systematic *>::iterator it = systMap.begin(); it != systMap.end(); ++it)
    lh.AddSystematic(it->second, systGroup[it->first]);
  // Add our pdfs
  lh.AddPdfs(pdfs, pdfGroups, genRates, norm_fitting_statuses);

  // And constraints
  std::vector<std::string> corrPairs;
  for (ParameterDict::iterator it = constrCorrs.begin(); it != constrCorrs.end(); ++it)
  {
    lh.SetConstraint(it->first, constrMeans.at(it->first), constrSigmas.at(it->first), constrCorrParName.at(it->first),
                     constrMeans.at(constrCorrParName.at(it->first)), constrSigmas.at(constrCorrParName.at(it->first)), it->second);
    corrPairs.push_back(constrCorrParName.at(it->first));
  }
  for (ParameterDict::iterator it = constrMeans.begin(); it != constrMeans.end(); ++it)
  {
    // Only add single parameter constraint if correlation hasn't already been applied
    if (constrCorrs.find(it->first) == constrCorrs.end() && std::find(corrPairs.begin(), corrPairs.end(), it->first) == corrPairs.end())
      lh.SetConstraint(it->first, it->second, constrSigmas.at(it->first));
  }
  for (ParameterDict::iterator it = constrRatioMeans.begin(); it != constrRatioMeans.end(); ++it)
    lh.SetConstraint(it->first, constrRatioParName.at(it->first), it->second, constrRatioSigmas.at(it->first));

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

  // Calculate the nominal LLH
  double nomllh = lh.Evaluate();

  if (isFakeData)
  {
    for (ParameterDict::iterator it = mins.begin(); it != mins.end(); ++it)
      parameterValues[it->first] = fdValues[it->first];

    lh.SetParameters(parameterValues);
    nomllh = lh.Evaluate();
  }

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

    // And a bit of jiggery pokery here to guarantee that the nominal value is one of the scan points
    double width = (maxs[name] - mins[name]) / (npoints);
    int numStepsBelowNom = floor((noms[name] - mins[name]) / width);
    int numStepsAboveNom = floor((maxs[name] - noms[name]) / width);

    double min = noms[name] - numStepsBelowNom * width;
    double max = noms[name] + numStepsAboveNom * width;

    // Make histogram for this parameter
    TString htitle = Form("%s, Asimov Rate: %f", name.c_str(), nom);
    TH1D *hScan = new TH1D((name + "_full").c_str(), (name + "_full").c_str(), npoints, (min - (width / 2)) / nom, (max + (width / 2)) / nom);
    hScan->SetTitle(std::string(htitle + ";" + name + " (rel. to Asimov); -(ln L_{full})").c_str());

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
      // Set bin contents
      hScan->SetBinContent(i + 1, llh - nomllh);

      // Return to nominal value
      parameterValues[name] = tempval;
    }
    // Write Histos
    hScan->Write();
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
