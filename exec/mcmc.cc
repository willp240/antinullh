// Sample the floating-oscillation model with Metropolis MCMC
// Reactor oscillation varies through the configured shape systematic
// Antinu headers
#include <AnalysisConfig.hh>
#include <DatasetSetup.hh>
#include <OutputDirs.hh>
#include <FitSummary.hh>
#include <FitParameters.hh>
#include <FitConstraints.hh>
#include <CombineLLHs.hh>
#include <BuildDatasets.hh>
#include <BuildComponent.hh>
#include <BuildLikelihood.hh>
#include <WriteDists.hh>
#include <RateReporting.hh>

// OXO headers
#include <StatisticSum.h>
#include <MCMC.h>
#include <MetropolisSampler.h>
#include <Rand.h>
#include <IO.h>

// ROOT headers
#include <TFile.h>
#include <TStopwatch.h>
#include <fstream>

using namespace antinufit;

void mcmc(const std::string &fitConfigFile_,
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
  const std::string projDir1D = dirs.outDir + "/1dlhproj";
  const std::string projDir2D = dirs.outDir + "/2dlhproj";
  if (fitConfig.GetSaveOutputs())
  {
    EnsureDirs(dirs);
    EnsureDir(projDir1D);
    EnsureDir(projDir2D);
  }

  // Build systematics
  SystsByDataset syst = BuildSysts(config.systematics,
                                    setup.pdfConfig,
                                    fitParameters.noms,
                                    setup.dsPDFMap,
                                    indexDistance);

  // Read nominal oscillation values and identify the theta12 parameter convention
  OscParams osc = BuildOscParams(fitParameters.noms);

  // Remove fit parameters not used by any dataset
  DatasetParameterMap datasetPars = PruneUnusedParameters(fitParameters, setup, config.systematics, osc);

  PrintFitSummary(fitParameters, datasetPars, osc, setup.pdfConfig);

  // Build component PDFs, the Asimov distribution, and the data distribution for each dataset
  PreparedDatasetMap datasets = BuildDatasets(fitConfig, fitParameters, setup, syst, config.systematics, osc, indexDistance, datasetPars);
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

  // Sampling settings remain specific to MCMC
  const int nsteps = fitConfig.GetIterations();
  const int burnin = fitConfig.GetBurnIn();
  const double sigmaScale = fitConfig.GetSigmaScale();
  const std::string &outDir = dirs.outDir;
  const std::string &postfitDistDir = dirs.postfitDistDir;
  ParameterDict &mins = fitParameters.mins;
  ParameterDict &maxs = fitParameters.maxs;
  ParameterDict &noms = fitParameters.noms;
  ParameterDict sigmas = fitParameters.sigmas;
  const ParameterDict &nbins = fitParameters.nbins;
  const ParameterDict &constrSigmas = fitParameters.constraintSigmas;
  std::map<std::string, bool> &fixedPars = fitParameters.fixedPars;

  // Starting point:
  //  - a fixed parameter stays nominal (the sampler never moves it);
  //  - the oscillation parameters, and any free parameter with no Gaussian
  //    prior on itself, start uniformly in [min, max];
  //  - every other free parameter starts at nominal + Gaussian(0, prior sigma), redrawn until inside [min, max]
  const std::set<std::string> oscParNames = {"deltam21", osc.theta12name};

  ParameterDict startVals = allParVals;
  std::cout << "Starting point:" << std::endl;
  for (ParameterDict::iterator parIt = startVals.begin(); parIt != startVals.end(); ++parIt)
  {
    const std::string &name = parIt->first;
    std::string rule;
    if (fixedPars.count(name) && fixedPars.at(name))
      rule = "fixed: nominal";
    else if (oscParNames.count(name) || !constrSigmas.count(name))
    {
      parIt->second = Rand::UniformRange(mins.at(name), maxs.at(name));
      rule = oscParNames.count(name) ? "oscillation: uniform in [min, max]" : "free, no Gaussian prior: uniform in [min, max]";
    }
    else
    {
      for (int iDraw = 0; iDraw < 10000; iDraw++)
      {
        const double draw = noms.at(name) + Rand::Gaus(0, constrSigmas.at(name));
        if (draw >= mins.at(name) && draw <= maxs.at(name))
        {
          parIt->second = draw;
          break;
        }
      }
      rule = "nominal + Gaussian(0, prior sigma)";
    }
    std::cout << "  " << name << " = " << parIt->second << "  [" << rule << "]" << std::endl;
  }

  fullLLH.SetParameters(startVals);
  TStopwatch timer;
  timer.Start(true);
  double startLLH = fullLLH.Evaluate();
  timer.Stop();
  std::cout << "NLLH at the starting point: " << startLLH << ", eval time: " << timer.RealTime() << " s" << std::endl;

  // Scale the configured sigmas to set the width of each Gaussian MCMC proposal
  MetropolisSampler sampler;
  for (ParameterDict::iterator it = sigmas.begin(); it != sigmas.end(); ++it)
    sigmas[it->first] *= sigmaScale;
  sampler.SetSigmas(sigmas);
  // Keep fixed parameters at their current values when proposing a new step
  for (std::map<std::string, bool>::iterator fixIt = fixedPars.begin(); fixIt != fixedPars.end(); ++fixIt)
    if (fixIt->second)
      sampler.Fix(fixIt->first);

  MCMC mh(sampler);
  mh.SetSaveChain(fitConfig.GetSaveOutputs());
  mh.SetMaxIter(nsteps);
  mh.SetBurnIn(burnin);
  mh.SetMinima(mins);
  mh.SetMaxima(maxs);
  mh.SetInitialTrial(startVals);

  // We're going to minmise not maximise -log(lh) and the
  // mc chain needs to know that the test stat is logged
  // otherwise it will give us the distribution of the log(lh) not the lh
  mh.SetTestStatLogged(true);
  mh.SetFlipSign(true);

  // Use configured bounds and bin counts for posterior histograms and projections
  AxisCollection lhAxes;
  for (ParameterDict::iterator it = mins.begin(); it != mins.end(); ++it)
    lhAxes.AddAxis(BinAxis(it->first, it->second, maxs.at(it->first), nbins.at(it->first)));
  mh.SetHistogramAxes(lhAxes);

  TStopwatch mcmcTimer;
  mcmcTimer.Start(true);
  const FitResult &res = mh.Optimise(&fullLLH);
  mcmcTimer.Stop();
  std::cout << "MCMC time: " << mcmcTimer.RealTime() << " s for " << nsteps << " steps (" << mcmcTimer.RealTime() / std::max(nsteps, 1) << " s/step)" << std::endl;

  ParameterDict bestFit = res.GetBestFit();
  fullLLH.SetParameters(bestFit);
  const double finalLLH = fullLLH.Evaluate();
  const bool validFit = res.GetValid();
  std::cout << "LLH: " << finalLLH << "\nFit Valid: " << validFit << std::endl;

  // Now save the results
  if (fitConfig.GetSaveOutputs())
  {
    res.SaveAs(outDir + "/fit_result.txt");
    std::ofstream file(outDir + "/fit_result.txt", std::ios::app);
    file << "\nLLH: " << finalLLH << "\n";
    file << "\nFit Valid: " << validFit << std::endl;
    std::cout << "Saved fit result to " << outDir + "/fit_result.txt" << std::endl;
  }

  MCMCSamples samples = mh.GetSamples();

  // Save the output tree
  if (fitConfig.GetSaveOutputs())
  {
    // Output filename is aboveDirectory_resultDirectory.root, inside outDir
    std::string chainFileName = outDir;
    std::string tempString1 = outDir;
    std::string tempString2 = outDir;
    size_t last_slash_idx = tempString1.find_last_of("/");
    tempString2.erase(0, last_slash_idx + 1);
    if (std::string::npos != last_slash_idx)
    {
      tempString1.erase(last_slash_idx, std::string::npos);
      last_slash_idx = tempString1.find_last_of("/");
      if (std::string::npos != last_slash_idx)
        tempString1.erase(0, last_slash_idx);
    }
    chainFileName = outDir + tempString1 + "_" + tempString2 + ".root";
    TFile *f = new TFile(chainFileName.c_str(), "recreate");
    TTree *outchain = samples.GetChain();
    outchain->Write();
    delete f;
  }

  // Save the histograms, and print the posterior mean and standard deviation of each free parameter
  typedef std::map<std::string, Histogram> HistMap;
  const HistMap &proj1D = samples.Get1DProjections();
  const HistMap &proj2D = samples.Get2DProjections();

  if (fitConfig.GetSaveOutputs())
    std::cout << "Saving LH projections to \n\t" << projDir1D << "\n\t" << projDir2D << std::endl;
  std::cout << "Posterior (after burn-in, from the 1D projections): mean +- std" << std::endl;
  for (HistMap::const_iterator it = proj1D.begin(); it != proj1D.end(); ++it)
  {
    if (fitConfig.GetSaveOutputs())
      IO::SaveHistogram(it->second, projDir1D + "/" + it->first + ".root");
    if (fixedPars[it->first])
      continue;
    const BinAxis &axis = it->second.GetAxes().GetAxis(0);
    double sum = 0, sumX = 0, sumXX = 0;
    for (size_t iBin = 0; iBin < axis.GetNBins(); iBin++)
    {
      const double w = it->second.GetBinContent(iBin);
      sum += w;
      sumX += w * axis.GetBinCentre(iBin);
      sumXX += w * axis.GetBinCentre(iBin) * axis.GetBinCentre(iBin);
    }
    if (sum > 0)
      std::cout << "  " << it->first << ": " << sumX / sum << " +- " << sqrt(std::max(0.0, sumXX / sum - sumX * sumX / sum / sum))
                << "   (nominal " << noms[it->first] << ")" << std::endl;
  }
  if (fitConfig.GetSaveOutputs())
    for (HistMap::const_iterator it = proj2D.begin(); it != proj2D.end(); ++it)
      IO::SaveHistogram(it->second, projDir2D + "/" + it->first + ".root");

  // Rebuild each dataset's predicted distribution using the fitted values
  ComponentRates bestFitRates;
  for (PreparedDatasetMap::const_iterator dataset = datasets.begin();
       dataset != datasets.end(); ++dataset)
  {
    const std::string &dsName = dataset->first;
    BinnedED total("postfit dist", setup.systAxes);
    total.SetObservables(setup.dataObs);
    // Reactor rate factors are already included in bestFit, so use 1.0 here
    std::vector<BinnedED> components = BuildComponents(
        dataset->second.pdfs, dataset->second.pdfGroups,
        syst, config.systematics, dsName, bestFit,
        std::vector<double>(dataset->second.pdfs.size(), 1.0), total, setup.dataObs);
    // Count events outside the likelihood buffer bins for the rate summary
    bestFitRates[dsName] = GetComponentRates(
        dataset->second.pdfs, components,
        setup.pdfConfig, setup.dataObs);
    if (fitConfig.GetSaveOutputs())
      WritePostfitDists(dsName, dataset->second, setup, components, total, dirs);
  }
  // Compare these fitted rates with the nominal rates saved during dataset building
  ComponentRates nominalRates;
  for (PreparedDatasetMap::const_iterator dataset = datasets.begin(); dataset != datasets.end(); ++dataset)
    nominalRates[dataset->first] = dataset->second.nominalNonBufferRates;
  PrintNominalBestFitSummary(fitParameters, bestFit, nominalRates, bestFitRates);

  if (fitConfig.GetSaveOutputs())
  {
    // Save autocorrelations
    std::ofstream cofs((outDir + "/auto_correlations.txt").c_str());
    std::vector<double> autocors = samples.GetAutoCorrelations();
    for (size_t i = 0; i < autocors.size(); i++)
      cofs << i << "\t" << autocors.at(i) << "\n";
    cofs.close();
  }

  std::cout << "Fit complete" << std::endl;
}

int main(int argc, char *argv[])
{
  if (argc != 5)
  {
    std::cout << "\nUsage: mcmc <fit_config_file> <eve_config_file> <pdf_config_file> <syst_config_file>" << std::endl;
    return 1;
  }

  std::string fitConfigFile(argv[1]);
  std::string eveConfigFile(argv[2]);
  std::string pdfConfigPath(argv[3]);
  std::string systConfigFile(argv[4]);

  mcmc(fitConfigFile, eveConfigFile, pdfConfigPath, systConfigFile);

  return 0;
}
