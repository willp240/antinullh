#include <string>
#include <FitConfigLoader.hh>
#include <DistConfigLoader.hh>
#include <DistConfig.hh>
#include <DistBuilder.hh>
#include <FitConfigLoader.hh>
#include <FitConfig.hh>
#include <EventConfigLoader.hh>
#include <EventConfig.hh>
#include <SystConfigLoader.hh>
#include <SystConfig.hh>
#include <SystFactory.hh>
#include <fstream>
#include <ROOTNtuple.h>
#include <BinnedNLLH.h>
#include <sys/stat.h>
#include <Rand.h>
#include <AxisCollection.h>
#include <IO.h>
#include <TH1D.h>
#include <Scale.h>

using namespace antinufit;

void llh_scan(const std::string &mcmcConfigFile_,
              const std::string &evConfigFile_,
              const std::string &pdfConfigFile_,
              const std::string &systConfigFile_)
{
  Rand::SetSeed(0);

  // Load up the fit configuration information
  FitConfig mcConfig;
  FitConfigLoader mcLoader(mcmcConfigFile_);
  mcConfig = mcLoader.LoadActive();
  bool isAsimov = mcConfig.GetAsimov();
  double livetime = mcConfig.GetLivetime();
  bool beestonBarlowFlag = mcConfig.GetBeestonBarlow();
  std::string outDir = mcConfig.GetOutDir();

  struct stat st = {0};
  if (stat(outDir.c_str(), &st) == -1)
    mkdir(outDir.c_str(), 0700);

  // Load up all the event types we want to contribute
  typedef std::map<std::string, EventConfig> EvMap;
  EventConfigLoader loader(evConfigFile_);
  EvMap toGet = loader.LoadActive();

  // Load up the PDF information (skeleton axis details, rather than the distributions themselves)
  DistConfigLoader pdfLoader(pdfConfigFile_);
  DistConfig pdfConfig = pdfLoader.Load();
  std::string pdfDir = pdfConfig.GetPDFDir();
  std::vector<std::string> dataObs = pdfConfig.GetBranchNames();
  ObsSet dataObsSet(dataObs);

  // Load up the systematics
  SystConfigLoader systLoader(systConfigFile_);
  SystConfig systConfig = systLoader.LoadActive();
  AxisCollection systAxes = DistBuilder::BuildAxes(pdfConfig);
  ParameterDict systNom = systConfig.GetNominal();
  ParameterDict systMaxima = systConfig.GetMaxima();
  ParameterDict systMinima = systConfig.GetMinima();
  ParameterDict systMass = systConfig.GetMass();
  ParameterDict systNbins = systConfig.GetNBins();
  std::map<std::string, std::string> systGroup = systConfig.GetGroup();
  std::map<std::string, std::string> systType = systConfig.GetType();
  std::map<std::string, std::string> systObs = systConfig.GetObs();

  // Loop over systematics and declare each one
  std::map<std::string, Systematic *> systMap;
  for (std::map<std::string, std::string>::iterator it = systType.begin(); it != systType.end(); ++it)
  {
    std::vector<std::string> obs;
    obs.push_back(systObs[it->first]);
    ObsSet obsSet(obs);

    Systematic *syst = SystFactory::New(it->first, systType[it->first], systNom[it->first]);
    syst->SetAxes(systAxes);
    // The "dimensions" the systematic applies too
    syst->SetTransformationObs(obsSet);
    // All the "dimensions" of the dataset
    syst->SetDistributionObs(dataObsSet);
    syst->Construct();
    systMap[it->first] = syst;
  }

  // Create the individual PDFs and Asimov components (could these be maps not vectors?)
  std::vector<BinnedED> pdfs;
  std::vector<BinnedED> indivAsmvDists;
  ParameterDict asimovRates;
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
      // WP check the pruned here
      std::cout << it->second.GetPrunedPath() << std::endl;
      dataSet = new ROOTNtuple(it->second.GetPrunedPath(), "pruned");
    }
    catch (const IOError &e_)
    {
      std::cout << "Warning: skipping " << it->first << " couldn't open files:\n\t" << e_.what() << std::endl;
      continue;
    }

    // Build distribution of those events
    BinnedED dist;
    dist = DistBuilder::Build(it->first, pdfConfig, dataSet);

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

    // Let's save the PDF as a histo
    IO::SaveHistogram(dist.GetHistogram(), pdfDir + "/" + it->first + ".root", dist.GetName());

    // Now scale the Asimov component by expected rate
    double rate = it->second.GetRate();
    dist.Scale(livetime * rate);
    asimov.Add(dist);
    indivAsmvDists.push_back(dist);
    asimovRates[it->first] = dist.Integral();
  }

  // And save combined histogram (final asimov dataset)
  IO::SaveHistogram(asimov.GetHistogram(), outDir + "/asimov.root", "asimov");

  // Now let's load up the data
  BinnedED dataDist;
  if (!isAsimov)
  {

    std::string dataPath = mcConfig.GetDatafile();

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
  lh.SetBufferAsOverflow(true);
  // lh.SetBuffer("energy", 1, 1);
  // Add our PDFs and data
  lh.AddPdfs(pdfs, pdfGroups, genRates);
  lh.SetDataDist(dataDist);
  // Set whether or not to use Beeston Barlow
  lh.SetBarlowBeeston(beestonBarlowFlag);
  // Add the systematics and any prior constraints
  for (std::map<std::string, Systematic *>::iterator it = systMap.begin(); it != systMap.end(); ++it)
    lh.AddSystematic(it->second, systGroup[it->first]);
  ParameterDict constrMeans = mcConfig.GetConstrMeans();
  ParameterDict constrSigmas = mcConfig.GetConstrSigmas();
  for (ParameterDict::iterator it = constrMeans.begin(); it != constrMeans.end(); ++it)
    lh.SetConstraint(it->first, it->second, constrSigmas.at(it->first));

  // Now add max and min values for each parameter
  ParameterDict mins = mcConfig.GetMinima();
  ParameterDict maxs = mcConfig.GetMaxima();
  for (ParameterDict::iterator it = systMass.begin(); it != systMass.end(); ++it)
  {
    mins[it->first] = systMinima[it->first];
    maxs[it->first] = systMaxima[it->first];
    asimovRates[it->first] = systNom[it->first];
  }

  // And finally bring it all together
  lh.RegisterFitComponents();

  // Now onto the LLH Scan. First set the number of points in the scan
  int npoints = 150;
  int countwidth = double(npoints) / double(5);

  // Initialise to nominal values
  ParameterDict parameterValues;
  for (ParameterDict::iterator it = mins.begin(); it != mins.end(); ++it)
    parameterValues[it->first] = asimovRates[it->first];

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
    double nom = asimovRates[name];
    if (nom == 0)
      nom = 1;
    double min = mins[name];
    double max = maxs[name];
    std::cout << "Scanning for " << name << " from " << min << " to " << max << ". Nom " << nom << std::endl;
    
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
  if (argc != 5)
  {
    std::cout << "\nUsage: llh_scan <mcmc_config_file> <eve_config_file> <pdf_config_file> <syst_config_file>" << std::endl;
    return 1;
  }

  std::string fitConfigFile(argv[1]);
  std::string eveConfigFile(argv[2]);
  std::string pdfConfigPath(argv[3]);
  std::string systConfigFile(argv[4]);

  llh_scan(fitConfigFile, eveConfigFile, pdfConfigPath, systConfigFile);

  return 0;
}
