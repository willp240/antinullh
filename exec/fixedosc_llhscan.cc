// Antinu headers
#include <FitConfigLoader.hh>
#include <OscGridConfigLoader.hh>
#include <DatasetSetup.hh>
#include <OutputDirs.hh>
#include <FitSummary.hh>
#include <ParamUtils.hh>
#include <CombineLLHs.hh>
#include <BuildDatasets.hh>

// OXO headers
#include <IO.h>
#include <Rand.h>

// ROOT headers
#include <TH1D.h>

#include <cmath>

using namespace antinufit;

namespace
{
  double EvaluateScanDeltaLLH(const BuildResult &datasets,
                              const SystSetup &syst,
                              const FitInputs &in,
                              size_t scanIndex,
                              double nomllh)
  {
    std::vector<BinnedNLLH> oscTestStats = BuildScanPointLLHs(datasets, syst, in, scanIndex);
    CombinedLLH oscCombined = CombineLLHs(oscTestStats, datasets.initialByDataset);
    StatisticSum &fullOscLLH = oscCombined.llh;
    ApplyConstraints(fullOscLLH, in.constraints, in.p.noms);
    fullOscLLH.RegisterFitComponents();
    return fullOscLLH.Evaluate() - nomllh;
  }
} // namespace

void fixedosc_llhscan(const std::string &fitConfigFile_,
                      const std::string &evConfigFile_,
                      const std::string &pdfConfigFile_,
                      const std::string &systConfigFile_,
                      const std::string &oscGridConfigFile_)
{
  Rand::SetSeed(0);

  // Load up the fit configuration information
  FitConfigLoader fitLoader(fitConfigFile_);
  FitConfig fitConfig = fitLoader.LoadActive();
  FitInputs in = MakeFitInputs(fitConfig);

  // Make sure output dirs exist
  OutputDirs dirs = MakeOutputDirs(in.run.outDir);
  if (in.run.saveOutputs)
    EnsureDirs(dirs);

  // Load in pdf axis and data path info
  DatasetSetup setup = LoadDatasetSetup(evConfigFile_, pdfConfigFile_);

  // Load up the oscillation probability grids
  OscGridConfigLoader oscGridLoader(oscGridConfigFile_);
  OscGridConfig oscGridConfig = oscGridLoader.Load();
  // Read the reactor distance info
  std::string reactorjson = oscGridConfig.GetReactorsJsonFile();
  std::unordered_map<int, double> indexDistance = LoadIndexDistanceMap(reactorjson);

  // Build systematics store
  SystSetup syst = BuildSystematics(systConfigFile_,
                                    setup.pdfConfig,
                                    in.p.noms,
                                    setup.dsPDFMap,
                                    indexDistance);

  // Build osc pars
  OscParams osc = BuildOscParams(in.p.noms);

  // Define the number of points
  int npoints = 150;
  int countwidth = double(npoints) / double(5);

  // Prune any parameters that aren't in fit config
  DatasetParsInfo dsParsInfo = PruneUnusedParameters(in, setup, syst, osc);

  PrintFitSummary(in, dsParsInfo.datasetPars, osc);

  // Legacy aliases used by scan code below
  const std::string outDir = in.run.outDir;
  const ParameterDict preBuildMins = in.p.mins;
  const ParameterDict preBuildMaxs = in.p.maxs;
  const ParameterDict preBuildNoms = in.p.noms;
  std::map<std::string, std::string> labelName = fitConfig.GetTexLabels();

  const std::string theta12name = osc.theta12name;
  const double deltam21_nom = osc.deltam21;
  const double deltam21_min = preBuildMins.at("deltam21");
  const double deltam21_max = preBuildMaxs.at("deltam21");
  const double theta12_nom = preBuildNoms.at(theta12name);
  const double theta12_min = preBuildMins.at(theta12name);
  const double theta12_max = preBuildMaxs.at(theta12name);

  // Match historical scan point placement: force nominal value to land on a scan point.
  const double dmWidth = (deltam21_max - deltam21_min) / static_cast<double>(npoints);
  const int dmStepsBelowNom = static_cast<int>(floor((deltam21_nom - deltam21_min) / dmWidth));
  const double dmScanMin = deltam21_nom - dmStepsBelowNom * dmWidth;

  const double thWidth = (theta12_max - theta12_min) / static_cast<double>(npoints);
  const int thStepsBelowNom = static_cast<int>(floor((theta12_nom - theta12_min) / thWidth));
  const double thScanMin = theta12_nom - thStepsBelowNom * thWidth;

  BuildDatasetsConfig buildCfg;
  buildCfg.component.buildScanVariants = true;

  // First npoints: deltam21 scan at nominal theta12
  for (int i = 0; i < npoints; ++i)
  {
    const double dm = dmScanMin + static_cast<double>(i) * dmWidth;

    OscillationPoint p;
    p.deltam21 = dm;
    p.theta12SinSq = Theta12ToSinSq(theta12name, theta12_nom);
    buildCfg.component.scanPoints.push_back(p);
  }

  // Second npoints: theta12 scan at nominal deltam21
  for (int i = 0; i < npoints; ++i)
  {
    const double th = thScanMin + static_cast<double>(i) * thWidth;

    OscillationPoint p;
    p.deltam21 = deltam21_nom;
    p.theta12SinSq = Theta12ToSinSq(theta12name, th);
    buildCfg.component.scanPoints.push_back(p);
  }

  BuildResult datasets = BuildDatasets(
      in, dirs, setup, syst, osc, indexDistance, dsParsInfo, buildCfg);

  // Use post-build values to match legacy behaviour (reactor scaling mutates these).
  ParameterDict mins = in.p.mins;
  ParameterDict maxs = in.p.maxs;
  ParameterDict noms = in.p.noms;

  CombinedLLH combined = CombineLLHs(datasets.llhs, datasets.initialByDataset);
  StatisticSum &fullLLH = combined.llh;
  ParameterDict allParVals = combined.initialValues;

  ApplyConstraints(fullLLH, in.constraints, in.p.noms);
  fullLLH.RegisterFitComponents();

  // Calculate the nominal LLH
  double nomllh = fullLLH.Evaluate();

  // Setup outfile
  std::string outFileName = outDir + "/llh_scan.root";
  TFile *outFile = new TFile(outFileName.c_str(), "recreate");

  // Loop over (non-oscillation) parameters
  for (ParameterDict::iterator parIt = mins.begin(); parIt != mins.end(); ++parIt)
  {
    // Get param name
    std::string name = parIt->first;

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
    TString htitle = Form("%s, Asimov Rate: %f", labelName[name].c_str(), nom);
    TH1D *hScan = new TH1D((name + "_full").c_str(), (labelName[name] + "_full").c_str(), npoints, (min - (width / 2)) / nom, (max + (width / 2)) / nom);
    hScan->SetTitle(std::string(htitle + ";" + labelName[name] + " (rel. to Asimov); -(ln L_{full})").c_str());

    std::cout << "Scanning for " << name << std::endl;
    // Now loop from min to max in npoint steps
    for (int i = 0; i < npoints; i++)
    {

      if (i % countwidth == 0)
        std::cout << i << "/" << npoints << " (" << double(i) / double(npoints) * 100 << "%)" << std::endl;

      // Set Parameters
      double parval = hScan->GetBinCenter(i + 1) * nom;
      double tempval = allParVals[name];
      allParVals[name] = parval;
      fullLLH.SetParameters(allParVals);

      // Evaluate LLH (later do sample and penalty)
      double llh = fullLLH.Evaluate();
      // Set bin contents
      hScan->SetBinContent(i + 1, llh - nomllh);
      std::cout << name << " " << parval << " " << llh - nomllh << std::endl;

      // Return to nominal value
      allParVals[name] = tempval;
      fullLLH.SetParameters(allParVals);
    }
    // Write Histos
    hScan->Write();
  }

  // Now we do the same for the oscillation parameters

  // And a bit of jiggery pokery here to guarantee that the nominal value is one of the scan points
  double width = (deltam21_max - deltam21_min) / (npoints);
  int numStepsBelowNom = floor((deltam21_nom - deltam21_min) / width);
  int numStepsAboveNom = floor((deltam21_max - deltam21_nom) / width);
  double min = deltam21_nom - numStepsBelowNom * width;
  double max = deltam21_nom + numStepsAboveNom * width;

  TString htitle = Form("%s, Nom. Value: %f", labelName["deltam21"].c_str(), deltam21_nom);
  TH1D *hDeltam = new TH1D("deltam21_full", (labelName["deltam21"] + "_full").c_str(), npoints, (min - (width / 2)), (max + (width / 2)));
  hDeltam->SetTitle(std::string(htitle + "; " + labelName["deltam21"] + " (eV^{2}); -(ln L_{full})").c_str());

  std::cout << "Scanning for deltam21" << std::endl;
  for (int iDeltaM = 0; iDeltaM < npoints; iDeltaM++)
  {
    if (iDeltaM % countwidth == 0)
      std::cout << iDeltaM << "/" << npoints << " (" << double(iDeltaM) / double(npoints) * 100 << "%)" << std::endl;

    const double deltaLLH = EvaluateScanDeltaLLH(
        datasets, syst, in, static_cast<size_t>(iDeltaM), nomllh);
    hDeltam->SetBinContent(iDeltaM + 1, deltaLLH);
  }
  hDeltam->Write();

  // And a bit of jiggery pokery here to guarantee that the nominal value is one of the scan points
  width = (theta12_max - theta12_min) / (npoints);
  numStepsBelowNom = floor((theta12_nom - theta12_min) / width);
  numStepsAboveNom = floor((theta12_max - theta12_nom) / width);
  min = theta12_nom - numStepsBelowNom * width;
  max = theta12_nom + numStepsAboveNom * width;

  // Repeat for theta
  htitle = Form("%s, Nom. Value: %f", labelName[theta12name].c_str(), theta12_nom);
  TH1D *hTheta12 = new TH1D((theta12name + "_full").c_str(), (labelName[theta12name] + "_nom_full").c_str(), npoints, (min - (width / 2)), (max + (width / 2)));
  hTheta12->SetTitle(std::string(htitle + "; " + labelName[theta12name] + " (^{o}); -(ln L_{full})").c_str());

  std::cout << "Scanning for " << theta12name << std::endl;
  for (int iTheta12 = 0; iTheta12 < npoints; iTheta12++)
  {
    const size_t scanIndex = static_cast<size_t>(iTheta12 + npoints);

    if (iTheta12 % countwidth == 0)
      std::cout << iTheta12 << "/" << npoints << " (" << double(iTheta12) / double(npoints) * 100 << "%)" << std::endl;

    const double deltaLLH = EvaluateScanDeltaLLH(
        datasets, syst, in, scanIndex, nomllh);
    hTheta12->SetBinContent(iTheta12 + 1, deltaLLH);
  }
  hTheta12->Write();

  // Close file
  outFile->Close();
  std::cout << "Wrote scan to " << outFile->GetName() << std::endl;

  delete outFile;
}

int main(int argc, char *argv[])
{
  if (argc != 6)
  {
    std::cout << "\nUsage: fixedosc_llhscan <fit_config_file> <eve_config_file> <pdf_config_file> <syst_config_file> <oscgrid_config_file>" << std::endl;
    return 1;
  }

  std::string fitConfigFile(argv[1]);
  std::string eveConfigFile(argv[2]);
  std::string pdfConfigPath(argv[3]);
  std::string systConfigFile(argv[4]);
  std::string oscgridConfigFile(argv[5]);

  fixedosc_llhscan(fitConfigFile, eveConfigFile, pdfConfigPath, systConfigFile, oscgridConfigFile);

  return 0;
}
