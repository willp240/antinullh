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

  // Load up the configuration data
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

  // Load up the PDFs
  DistConfigLoader pdfLoader(pdfConfigFile_);
  DistConfig pdfConfig = pdfLoader.Load();
  std::string pdfDir = pdfConfig.GetPDFDir();

  // Load up the systematics
  SystConfigLoader systLoader(systConfigFile_);
  SystConfig systConfig = systLoader.LoadActive();

  AxisCollection systAxes = DistBuilder::BuildAxes(pdfConfig);
  std::vector<std::string> dataObs = pdfConfig.GetBranchNames();
  ObsSet dataObsSet(dataObs);

  ParameterDict syst_nom = systConfig.GetNominal();
  ParameterDict syst_maxima = systConfig.GetMaxima();
  ParameterDict syst_minima = systConfig.GetMinima();
  ParameterDict syst_mass = systConfig.GetMass();
  ParameterDict syst_nbins = systConfig.GetNBins();
  std::map<std::string, std::string> syst_type = systConfig.GetType();
  std::map<std::string, std::string> syst_obs = systConfig.GetObs();

  std::vector<Systematic *> syst_vec;
  // Loop over systematics and declare each type. Must be a better way to do this but it will do for now
  for (std::map<std::string, std::string>::iterator it = syst_type.begin(); it != syst_type.end(); ++it)
  {
    std::vector<std::string> Obs;
    Obs.push_back(syst_obs[it->first]);
    ObsSet obsSet(Obs);

    double stddev_nom = 0;
    if (syst_type[it->first] == "conv")
      stddev_nom = syst_nom[it->first + "_stddevs"];
    Systematic *syst = SystFactory::New(it->first, syst_type[it->first], syst_nom[it->first], stddev_nom);
    syst->SetAxes(systAxes);
    // The "dimensions" the systematic applies too
    syst->SetTransformationObs(obsSet);
    // All the "dimensions" of the dataset
    syst->SetDistributionObs(dataObsSet);
    syst->Construct();
    syst_vec.push_back(syst);
  }

  // Create the individual PDFs and Asimov components (could these be maps not vectors)
  std::vector<BinnedED> pdfs;
  std::vector<BinnedED> indivAsmvDists;
  ParameterDict asimovRates;
  std::vector<int> genRates;

  // Create the empty full dist
  BinnedED asimov = BinnedED("asimov", systAxes);
  asimov.SetObservables(dataObs);

  // Build each of the PDFs, scale them to the correct size
  for (EvMap::iterator it = toGet.begin(); it != toGet.end(); ++it)
  {

    std::cout << "Building distribution for " << it->first << std::endl;
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

    genRates.push_back(dist.Integral());

    // Scale for PDF
    if (dist.Integral())
      dist.Normalise();

    pdfs.push_back(dist);

    // Apply nominal systematic variables
    for (int i_syst = 0; i_syst < syst_vec.size(); i_syst++)
    {
      double distInt = dist.Integral();
      dist = syst_vec.at(i_syst)->operator()(dist);
      dist.Scale(distInt);
    }

    // Let's save the PDF as a histo
    IO::SaveHistogram(dist.GetHistogram(), pdfDir + "/" + it->first + ".root", dist.GetName());

    // Now scale the Asimov component by expected rate
    double rate = it->second.GetRate();
    dist.Scale(livetime * rate);
    std::cout << pdfs.back().GetBinContent(0) << std::endl;
    asimov.Add(dist);
    indivAsmvDists.push_back(dist);
    asimovRates[it->first] = dist.Integral();
  }

  IO::SaveHistogram(asimov.GetHistogram(), outDir + "/asimov.root", "asimov");

  // If its a root tree then load it up
  BinnedED dataDist;

  if (!isAsimov)
  {

    std::string dataPath = mcConfig.GetDatafile();

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
  lh.SetBufferAsOverflow(true);
  // lh.SetBuffer("energy", 1, 1);
  lh.AddPdfs(pdfs, genRates);
  lh.SetDataDist(dataDist);
  lh.SetBarlowBeeston(beestonBarlowFlag);
  for (int i_syst = 0; i_syst < syst_vec.size(); i_syst++)
    lh.AddSystematic(syst_vec.at(i_syst));

  ParameterDict constrMeans = mcConfig.GetConstrMeans();
  ParameterDict constrSigmas = mcConfig.GetConstrSigmas();
  ParameterDict mins = mcConfig.GetMinima();
  ParameterDict maxs = mcConfig.GetMaxima();

  for (ParameterDict::iterator it = constrMeans.begin(); it != constrMeans.end(); ++it)
    lh.SetConstraint(it->first, it->second, constrSigmas.at(it->first));

  for (ParameterDict::iterator it = syst_mass.begin(); it != syst_mass.end(); ++it)
  {
    mins[it->first] = syst_minima[it->first];
    maxs[it->first] = syst_maxima[it->first];
    asimovRates[it->first] = syst_nom[it->first];
  }

  lh.RegisterFitComponents();

  // number of points in scan
  int npoints = 150;
  int countwidth = 10;
  // double(npoints) / double(5);

  // Initialise to nominal values
  ParameterDict parameterValues;
  for (ParameterDict::iterator it = mins.begin(); it != mins.end(); ++it)
    parameterValues[it->first] = asimovRates[it->first];

  for (ParameterDict::iterator it = constrMeans.begin(); it != constrMeans.end(); ++it)
    parameterValues[it->first] = constrMeans[it->first];

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
    double nom = asimovRates[name];
    if (nom == 0)
      nom = 1;
    double min = mins[name];
    double max = maxs[name];
    std::cout << "Scanning for " << name << " from " << min << " to " << max << ". Nom " << nom << std::endl;
    // Make histos
    TString htitle = Form("%s, Asimov Rate: %f", name.c_str(), nom);
    TH1D *hScan = new TH1D((name + "_full").c_str(), (name + "_full").c_str(), npoints, min / nom, max / nom);
    hScan->SetTitle(std::string(htitle + ";" + name + " (rel. to Asimov); -(ln L_{full})").c_str());

    // Commented out bits here are in anticipation of splitting LLH calculation in oxo at some point. Could then see prior and sample contributions separately
    // TH1D *hScanSam = new TH1D((name+"_sam").c_str(), (name+"_sam").c_str(), npoints, min/nom, max/nom);
    // hScanSam->SetTitle(std::string(std::string("2LLH_sam, ") + name + ";" + name + "; -2(ln L_{sample})").c_str());
    // TH1D *hScanPen = new TH1D((name+"_pen").c_str(), (name+"_pen").c_str(), npoints, min/nom, max/nom);
    // hScanPen->SetTitle(std::string(std::string("2LLH_pen, ") + name + ";" + name + "; -2(ln L_{penalty})").c_str());

    // loop from min to max in steps of 150 (might want to do smaller range)
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
      double LLH = lh.Evaluate();

      // SetBinContents
      hScan->SetBinContent(i + 1, LLH);
      // hScanSam->SetBinContent(i+1, LLH);
      // hScanPen->SetBinContent(i+1, LLH);

      // return to asimov value
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
