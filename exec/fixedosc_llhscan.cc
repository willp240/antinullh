// Antinu headers
#include <AnalysisConfig.hh>
#include <DatasetSetup.hh>
#include <OutputDirs.hh>
#include <FitSummary.hh>
#include <FitParameters.hh>
#include <FitConstraints.hh>
#include <CombineLLHs.hh>
#include <BuildDatasets.hh>
#include <BuildDist.hh>
#include <LoadEvents.hh>
#include <BuildLikelihood.hh>
#include <WriteDists.hh>

// OXO headers
#include <IO.h>
#include <Rand.h>

// ROOT headers
#include <TH1D.h>
#include <TFile.h>

#include <algorithm>
#include <cmath>
#include <map>
#include <memory>
#include <stdexcept>
#include <utility>
#include <vector>

using namespace antinufit;

namespace
{
  struct OscillationPoint
  {
    double deltam21 = 0.0;
    double theta12SinSq = 0.0;
  };

  struct ScanComponent
  {
    BinnedED pdf;
    int generatedRate = 0;
    double reactorRatio = 1.0;
  };

  // dataset -> component name -> values for the current scan point
  using ScanPointData = std::map<std::string, std::map<std::string, ScanComponent>>;

  // Evaluate one oscillation point relative to the nominal LLH
  double EvaluateScanDeltaLLH(const PreparedDatasetMap &datasets,
                              const ScanPointData &pointData,
                              const SystsByDataset &syst,
                              const SystConfig &systConfig,
                              const FitConfig &fitConfig,
                              const FitParameters &fitParameters,
                              const PDFConfig &pdfConfig,
                              double nomllh)
  {
    std::vector<BinnedNLLH> pointLLHs;
    pointLLHs.reserve(datasets.size());
    std::map<std::string, double> factors;
    for (PreparedDatasetMap::const_iterator dataset = datasets.begin();
         dataset != datasets.end(); ++dataset)
    {
      const std::string &dsName = dataset->first;
      const PreparedDataset &prepared = dataset->second;
      // Keep non-reactor PDFs and replace only the reactor PDFs at this point
      std::vector<BinnedED> pdfs = prepared.pdfs;
      std::vector<int> generatedRates = prepared.genRates;

      ScanPointData::const_iterator scanDataset = pointData.find(dsName);
      if (scanDataset != pointData.end())
        for (std::map<std::string, ScanComponent>::const_iterator component =
                 scanDataset->second.begin(); component != scanDataset->second.end(); ++component)
        {
          const ScanComponent &scanComponent = component->second;
          size_t pdfIndex = 0;
          while (pdfIndex < pdfs.size() && pdfs[pdfIndex].GetName() != component->first)
            ++pdfIndex;
          if (pdfIndex == pdfs.size())
            throw std::runtime_error("Scan PDF " + component->first + " is not in dataset " + dsName);
          pdfs[pdfIndex] = scanComponent.pdf;
          generatedRates[pdfIndex] = scanComponent.generatedRate;

          // Adjust reactor rate parameters by the change from nominal survival
          const double baseRatio = prepared.reactorRatio.at(component->first);
          if (!(baseRatio > 0))
            throw std::runtime_error("Zero nominal reactor ratio for " + component->first);
          const double factor = scanComponent.reactorRatio / baseRatio;
          std::pair<std::map<std::string, double>::iterator, bool> inserted =
              factors.emplace(component->first, factor);
          if (!inserted.second && std::fabs(inserted.first->second - factor) > 1e-12)
            throw std::runtime_error("Shared reactor rate has different scan ratios: " + component->first);
        }

      // Build this dataset's likelihood with the selected reactor PDFs and counts
      BinnedNLLH &lh = pointLLHs.emplace_back();
      ConfigureDatasetLikelihood(lh, syst, systConfig, fitConfig, pdfConfig,
                                 dsName, prepared, pdfs, generatedRates);
    }
    StatisticSum fullOscLLH = CombineLLHs(pointLLHs);

    // Scale reactor rate values and Gaussian priors by the change from nominal
    ParameterDict values = GetInitialValues(datasets);
    ParameterDict nominal = fitParameters.noms;
    FitParameters scanParameters = fitParameters;
    for (std::map<std::string, double>::const_iterator entry = factors.begin();
         entry != factors.end(); ++entry)
    {
      values.at(entry->first) *= entry->second;
      nominal.at(entry->first) *= entry->second;
      if (scanParameters.constraintMeans.count(entry->first))
      {
        scanParameters.constraintMeans.at(entry->first) *= entry->second;
        scanParameters.constraintSigmas.at(entry->first) *= entry->second;
      }
    }
    // Evaluate the combined likelihood with the point-specific PDFs and priors
    ApplyConstraints(fullOscLLH, scanParameters, nominal);
    fullOscLLH.RegisterFitComponents();
    fullOscLLH.SetParameters(values);
    return fullOscLLH.Evaluate() - nomllh;
  }
} // namespace

// Scan each parameter while holding the others nominal and write delta-LLH histograms
void fixedosc_llhscan(const std::string &fitConfigFile_,
                      const std::string &evConfigFile_,
                      const std::string &pdfConfigFile_,
                      const std::string &systConfigFile_)
{
  Rand::SetSeed(0);

  // Load up the fit configuration information
  AnalysisConfig config = LoadAnalysisConfig(fitConfigFile_, evConfigFile_,
                                             pdfConfigFile_, systConfigFile_);
  const FitConfig &fitConfig = config.fit;
  const DatasetSetup &setup = config.datasets;
  const std::unordered_map<int, double> &indexDistance = config.reactorBaselines;
  // Copy parameter settings from the config so pruning and reactor scaling can modify them
  FitParameters fitParameters = MakeFitParameters(fitConfig);

  // Make sure output dirs exist
  OutputDirs dirs = MakeOutputDirs(fitConfig.GetOutDir());
  if (fitConfig.GetSaveOutputs())
    EnsureDirs(dirs);
  else
    EnsureDir(dirs.outDir); // The scan result is always written

  // Build systematics
  SystsByDataset syst = BuildSysts(config.systematics,
                                    setup.pdfConfig,
                                    fitParameters.noms,
                                    setup.dsPDFMap,
                                    indexDistance);

  // Read nominal oscillation values and identify the theta12 parameter convention
  OscParams osc = BuildOscParams(fitParameters.noms);

  // Define the number of points in the scan
  int npoints = fitConfig.GetLLHScanPoints();
  if (npoints <= 0)
  {
    std::cout << "WARNING: llh_scan_points <= 0 in fit config. Defaulting to 150." << std::endl;
    npoints = 150;
  }
  int countwidth = std::max(1, static_cast<int>(double(npoints) / double(5)));

  // Remove fit parameters not used by any dataset
  DatasetParameterMap datasetPars = PruneUnusedParameters(fitParameters, setup, config.systematics, osc);

  PrintFitSummary(fitParameters, datasetPars, osc, setup.pdfConfig);

  std::map<std::string, std::string> labelName = fitConfig.GetTexLabels();

  const std::string theta12name = osc.theta12name;
  const double deltam21_nom = osc.deltam21;
  // Use the configured oscillation bounds and nominal values to place scan points
  const double deltam21_min = fitParameters.mins.at("deltam21");
  const double deltam21_max = fitParameters.maxs.at("deltam21");
  const double theta12_nom = fitParameters.noms.at(theta12name);
  const double theta12_min = fitParameters.mins.at(theta12name);
  const double theta12_max = fitParameters.maxs.at(theta12name);

  // Start each scan grid so its nominal value is one of the scan points
  const double dmWidth = (deltam21_max - deltam21_min) / static_cast<double>(npoints);
  const int dmStepsBelowNom = static_cast<int>(floor((deltam21_nom - deltam21_min) / dmWidth));
  const double dmScanMin = deltam21_nom - dmStepsBelowNom * dmWidth;

  const double thWidth = (theta12_max - theta12_min) / static_cast<double>(npoints);
  const int thStepsBelowNom = static_cast<int>(floor((theta12_nom - theta12_min) / thWidth));
  const double thScanMin = theta12_nom - thStepsBelowNom * thWidth;

  std::vector<OscillationPoint> scanPoints;
  scanPoints.reserve(2 * static_cast<size_t>(npoints));

  // First npoints: deltam21 scan at nominal theta12
  for (int i = 0; i < npoints; ++i)
  {
    const double dm = dmScanMin + static_cast<double>(i) * dmWidth;

    OscillationPoint point;
    point.deltam21 = dm;
    point.theta12SinSq = Theta12ToSinSq(theta12name, theta12_nom);
    scanPoints.push_back(point);
  }

  // Second npoints: theta12 scan at nominal deltam21
  for (int i = 0; i < npoints; ++i)
  {
    const double th = thScanMin + static_cast<double>(i) * thWidth;

    OscillationPoint point;
    point.deltam21 = deltam21_nom;
    point.theta12SinSq = Theta12ToSinSq(theta12name, th);
    scanPoints.push_back(point);
  }

  // Build component PDFs, the Asimov distribution, and the data distribution for each dataset
  PreparedDatasetMap datasets = BuildDatasets(
      fitConfig, fitParameters, setup, syst, config.systematics, osc, indexDistance, datasetPars);
  // Use the calculated reactor ratios to scale fit rates, bounds, priors, and fake-data rates
  ApplyReactorRateScales(fitParameters, datasets);
  // Write input PDFs, Asimov distributions, and fake-data distributions if enabled
  if (fitConfig.GetSaveOutputs())
    WriteInputDists(datasets, setup, dirs);

  // Build one likelihood per dataset, then combine them
  std::vector<BinnedNLLH> datasetLLHs = BuildDatasetLikelihoods(datasets, syst, config.systematics, fitConfig, setup.pdfConfig);
  StatisticSum fullLLH = CombineLLHs(datasetLLHs);
  ParameterDict allParVals = GetInitialValues(datasets);

  // Finally apply all our constraints
  ApplyConstraints(fullLLH, fitParameters, fitParameters.noms);
  fullLLH.RegisterFitComponents();

  // Calculate the nominal LLH
  double nomllh = fullLLH.Evaluate();

  // Setup outfile
  std::string outFileName = dirs.outDir + "/llh_scan.root";
  TFile outFile(outFileName.c_str(), "RECREATE");
  if (outFile.IsZombie())
    throw std::runtime_error("Cannot create scan output " + outFileName);

  // Loop over (non-oscillation) parameters
  for (ParameterDict::const_iterator parIt = fitParameters.mins.begin();
       parIt != fitParameters.mins.end(); ++parIt)
  {
    // Get param name
    std::string name = parIt->first;
    if (name == "deltam21" || name == theta12name)
      continue;

    // Set scan range from this parameter's max and min values
    // We plot x axis as relative to the nominal value, so if this = 0 we set it to 1
    const double nominalValue = fitParameters.noms.at(name);
    double nom = nominalValue;
    if (nom == 0)
      nom = 1;
    else
      nom = std::abs(nom);

    // And a bit of jiggery pokery here to guarantee that the nominal value is one of the scan points
    const double lower = fitParameters.mins.at(name);
    const double upper = fitParameters.maxs.at(name);
    double width = (upper - lower) / (npoints);
    int numStepsBelowNom = floor((nominalValue - lower) / width);
    int numStepsAboveNom = floor((upper - nominalValue) / width);

    double min = nominalValue - numStepsBelowNom * width;
    double max = nominalValue + numStepsAboveNom * width;

    // Make histogram for this parameter
    TString htitle = Form("%s, Asimov Rate: %f", labelName[name].c_str(), nom);
    TH1D *hScan = new TH1D((name + "_full").c_str(), (labelName[name] + "_full").c_str(), npoints, (min - (width / 2)) / nom, (max + (width / 2)) / nom);
    hScan->SetTitle(std::string(htitle + ";" + labelName[name] + " (rel. to |Asimov|); -(ln L_{full})").c_str());

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

      // Evaluate LLH
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

  // Scan the oscillation parameters with a new reactor PDF at each point
  TString htitle = Form("%s, Nom. Value: %f", labelName["deltam21"].c_str(), deltam21_nom);
  TH1D *hDeltam = new TH1D("deltam21_full", (labelName["deltam21"] + "_full").c_str(),
                          npoints, dmScanMin - dmWidth / 2,
                          dmScanMin + (npoints - 0.5) * dmWidth);
  hDeltam->SetTitle(std::string(htitle + "; " + labelName["deltam21"] + " (eV^{2}); -(ln L_{full})").c_str());

  htitle = Form("%s, Nom. Value: %f", labelName[theta12name].c_str(), theta12_nom);
  TH1D *hTheta12 = new TH1D((theta12name + "_full").c_str(),
                            (labelName[theta12name] + "_nom_full").c_str(),
                            npoints, thScanMin - thWidth / 2,
                            thScanMin + (npoints - 0.5) * thWidth);
  hTheta12->SetTitle(std::string(htitle + "; " + labelName[theta12name] + " (^{o}); -(ln L_{full})").c_str());

  // Keep reactor event data open while visiting all scan points
  std::map<std::string, std::map<std::string, std::unique_ptr<DataSet>>> reactorEvents;
  for (size_t scanIndex = 0; scanIndex < scanPoints.size(); ++scanIndex)
  {
    const bool scanningDeltam = scanIndex < static_cast<size_t>(npoints);
    const size_t pointIndex = scanningDeltam ? scanIndex : scanIndex - npoints;
    if (pointIndex % countwidth == 0)
      std::cout << "Scanning " << (scanningDeltam ? "deltam21" : theta12name)
                << " " << pointIndex << "/" << npoints << std::endl;

    // Build the reactor PDFs and counts for this point before evaluating its LLH
    const OscillationPoint &point = scanPoints[scanIndex];
    ScanPointData pointData;
    for (DSMap::const_iterator dataset = setup.dsPDFMap.begin();
         dataset != setup.dsPDFMap.end(); ++dataset)
    {
      const std::string &dsName = dataset->first;
      const PreparedDataset &prepared = datasets.at(dsName);
      for (EvMap::const_iterator event = dataset->second.begin();
           event != dataset->second.end(); ++event)
      {
        if (!event->second.GetOscillated())
          continue;

        const std::string &pdfName = event->first;
        std::unique_ptr<DataSet> &events = reactorEvents[dsName][pdfName];
        if (!events)
        {
          events = LoadEvents(event->second.GetPrunedPath());
          if (!events)
            throw std::runtime_error("Cannot load MC for " + pdfName + " in " + dsName +
                                     ": " + event->second.GetPrunedPath());
        }

        size_t nominalIndex = 0;
        while (nominalIndex < prepared.pdfs.size() &&
               prepared.pdfs[nominalIndex].GetName() != pdfName)
          ++nominalIndex;
        if (nominalIndex == prepared.pdfs.size())
          throw std::runtime_error("Missing nominal PDF for " + pdfName + " in " + dsName);

        // Reuse the prepared PDF at the nominal point and rebuild it elsewhere
        ScanComponent component;
        component.reactorRatio = prepared.reactorRatio.at(pdfName);
        if (point.deltam21 == osc.deltam21 &&
            point.theta12SinSq == osc.theta12_sinsq())
        {
          component.pdf = prepared.pdfs[nominalIndex];
          component.generatedRate = prepared.genRates.at(nominalIndex);
        }
        else
        {
          component.pdf = BuildDist::BuildOscillatedDist(
              pdfName, event->second.GetNumDimensions(), setup.pdfConfig,
              events.get(), point.deltam21, point.theta12SinSq,
              indexDistance, component.reactorRatio);
          // Keep the generated count before normalising the likelihood PDF
          component.pdf.AddPadding();
          component.generatedRate = static_cast<int>(component.pdf.Integral());
          if (component.pdf.Integral())
            component.pdf.Normalise();
        }
        pointData[dsName][pdfName] = std::move(component);
      }
    }

    const double deltaLLH = EvaluateScanDeltaLLH(
        datasets, pointData, syst, config.systematics, fitConfig, fitParameters, setup.pdfConfig, nomllh);
    if (scanningDeltam)
      hDeltam->SetBinContent(pointIndex + 1, deltaLLH);
    else
      hTheta12->SetBinContent(pointIndex + 1, deltaLLH);
  }
  hDeltam->Write();
  hTheta12->Write();

  // Close file
  outFile.Close();
  std::cout << "Wrote scan to " << outFileName << std::endl;
}

// Read the four config paths and run the scan
int main(int argc, char *argv[])
{
  if (argc != 5)
  {
    std::cout << "\nUsage: fixedosc_llhscan <fit_config_file> <eve_config_file> <pdf_config_file> <syst_config_file>" << std::endl;
    return 1;
  }

  std::string fitConfigFile(argv[1]);
  std::string eveConfigFile(argv[2]);
  std::string pdfConfigPath(argv[3]);
  std::string systConfigFile(argv[4]);

  fixedosc_llhscan(fitConfigFile, eveConfigFile, pdfConfigPath, systConfigFile);

  return 0;
}
