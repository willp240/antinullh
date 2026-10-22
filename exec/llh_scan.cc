#include <AnalysisConfig.hh>
#include <DatasetSetup.hh>
#include <OutputDirs.hh>
#include <FitSummary.hh>
#include <FitParameters.hh>
#include <FitConstraints.hh>
#include <CombineLLHs.hh>
#include <BuildDatasets.hh>
#include <BuildLikelihood.hh>
#include <WriteDists.hh>

#include <Rand.h>
#include <TH1D.h>
#include <TFile.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>

using namespace antinufit;

void llh_scan(const std::string &fitConfigFile,
              const std::string &eventConfigFile,
              const std::string &pdfConfigFile,
              const std::string &systConfigFile)
{
  Rand::SetSeed(0);

  // Use the same floating-oscillation model and likelihood setup as mcmc
  AnalysisConfig config = LoadAnalysisConfig(fitConfigFile, eventConfigFile,
                                             pdfConfigFile, systConfigFile);
  const FitConfig &fitConfig = config.fit;
  const DatasetSetup &setup = config.datasets;
  const std::unordered_map<int, double> &baselines = config.reactorBaselines;
  // Copy parameter settings from the config so pruning and reactor scaling can modify them
  FitParameters fitParameters = MakeFitParameters(fitConfig);
  // Make sure output dirs exist
  const OutputDirs dirs = MakeOutputDirs(fitConfig.GetOutDir());
  if (fitConfig.GetSaveOutputs())
    EnsureDirs(dirs);
  else
    EnsureDir(dirs.outDir); // The scan itself is always the requested output

  // Build systematics
  SystsByDataset syst = BuildSysts(config.systematics, setup.pdfConfig,
                                    fitParameters.noms, setup.dsPDFMap, baselines);
  // Read nominal oscillation values and identify the theta12 parameter convention
  OscParams osc = BuildOscParams(fitParameters.noms);
  // Remove fit parameters not used by any dataset
  DatasetParameterMap datasetPars = PruneUnusedParameters(fitParameters, setup, config.systematics, osc);
  PrintFitSummary(fitParameters, datasetPars, osc, setup.pdfConfig);

  // Build component PDFs, the Asimov distribution, and the data distribution for each dataset
  PreparedDatasetMap datasets = BuildDatasets(fitConfig, fitParameters, setup, syst, config.systematics, osc,
                                       baselines, datasetPars);
  // Use the calculated reactor ratios to scale fit rates, bounds, priors, and fake-data rates
  ApplyReactorRateScales(fitParameters, datasets);

  // Write input PDFs, Asimov distributions, and fake-data distributions if enabled
  if (fitConfig.GetSaveOutputs())
    WriteInputDists(datasets, setup, dirs);

  // Build one likelihood per dataset, then combine them
  std::vector<BinnedNLLH> datasetLLHs = BuildDatasetLikelihoods(datasets, syst, config.systematics, fitConfig, setup.pdfConfig);
  StatisticSum llh = CombineLLHs(datasetLLHs);
  // Finally apply all our constraints
  ApplyConstraints(llh, fitParameters, fitParameters.noms);
  llh.RegisterFitComponents();

  // These are one-dimensional slices through the nominal point, not profiles
  const ParameterDict nominal = GetInitialValues(datasets);
  llh.SetParameters(nominal);
  const double nominalLLH = llh.Evaluate();

  int npoints = fitConfig.GetLLHScanPoints();
  if (npoints <= 0)
  {
    std::cout << "WARNING: llh_scan_points <= 0 in fit config. Defaulting to 150.\n";
    npoints = 150;
  }

  const std::string outputPath = dirs.outDir + "/llh_scan.root";
  TFile output(outputPath.c_str(), "RECREATE");
  if (output.IsZombie())
    throw std::runtime_error("Cannot create scan output " + outputPath);

  const std::map<std::string, std::string> labels = fitConfig.GetTexLabels();
  const int progressEvery = std::max(1, npoints / 5);
  for (ParameterDict::const_iterator entry = nominal.begin(); entry != nominal.end(); ++entry)
  {
    const std::string &name = entry->first;
    const double centre = entry->second;
    const double lower = fitParameters.mins.at(name);
    const double upper = fitParameters.maxs.at(name);
    const double step = (upper - lower) / npoints;
    if (!(step > 0))
      throw std::runtime_error("Invalid scan bounds for " + name);

    // Place the nominal point exactly at a bin centre, as in fixedosc_llhscan
    const int stepsBelow = static_cast<int>(std::floor((centre - lower) / step));
    const double first = centre - stepsBelow * step;
    const double scale = centre == 0.0 ? 1.0 : std::abs(centre);
    std::map<std::string, std::string>::const_iterator labelIt = labels.find(name);
    const std::string label = labelIt == labels.end() ? name : labelIt->second;
    TH1D scan((name + "_full").c_str(), label.c_str(), npoints,
              (first - step / 2) / scale,
              (first + (npoints - 0.5) * step) / scale);
    scan.SetTitle((label + ";" + label +
                   " (relative to |nominal|); -(ln L_{full})").c_str());

    std::cout << "Scanning " << name << '\n';
    ParameterDict trial = nominal;
    for (int i = 0; i < npoints; ++i)
    {
      if (i % progressEvery == 0)
        std::cout << i << "/" << npoints << '\n';
      trial[name] = first + i * step;
      llh.SetParameters(trial);
      scan.SetBinContent(i + 1, llh.Evaluate() - nominalLLH);
    }
    llh.SetParameters(nominal);
    scan.Write();
  }

  output.Close();
  std::cout << "Wrote scan to " << outputPath << '\n';
}

int main(int argc, char *argv[])
{
  if (argc != 5)
  {
    std::cerr << "Usage: llh_scan <fit_config_file> <event_config_file> "
                 "<pdf_config_file> <syst_config_file>\n";
    return 1;
  }
  llh_scan(argv[1], argv[2], argv[3], argv[4]);
  return 0;
}
