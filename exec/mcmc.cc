// Antinu headers
#include <FitConfigLoader.hh>
#include <DatasetSetup.hh>
#include <OutputDirs.hh>
#include <FitSummary.hh>
#include <ApplySystematics.hh>
#include <Projection.hh>
#include <ParamUtils.hh>
#include <CombineLLHs.hh>
#include <BuildDatasets.hh>
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

namespace
{
  bool IsFusedOsc(const std::string &type)
  {
    return type == "FusedOsc" || type.rfind("FusedOsc:", 0) == 0;
  }
}

void mcmc(const std::string &fitConfigFile_,
          const std::string &evConfigFile_,
          const std::string &pdfConfigFile_,
          const std::string &systConfigFile_)
{
  Rand::SetSeed(0);

  // Load up the fit configuration information
  FitConfigLoader fitLoader(fitConfigFile_);
  FitConfig fitConfig = fitLoader.LoadActive();
  FitInputs in = MakeFitInputs(fitConfig);

  // Make sure output dirs exist
  OutputDirs dirs = MakeOutputDirs(in.run.outDir);
  const std::string projDir1D = dirs.outDir + "/1dlhproj";
  const std::string projDir2D = dirs.outDir + "/2dlhproj";
  if (in.run.saveOutputs)
  {
    EnsureDirs(dirs);
    EnsureDir(projDir1D);
    EnsureDir(projDir2D);
  }

  // Load in pdf axis and data path info
  DatasetSetup setup = LoadDatasetSetup(evConfigFile_, pdfConfigFile_);

  // Reactor baselines are shared by fixed and floating oscillation calculations.
  const auto indexDistance = LoadIndexDistanceMap(fitConfig.GetReactorsJsonFile());

  SystSetup syst = BuildSystematics(systConfigFile_,
                                    setup.pdfConfig,
                                    in.p.noms,
                                    setup.dsPDFMap);

  // MCMC needs floating oscillation weights for every oscillated component.
  const bool hasFloatingOscillation = std::any_of(syst.type.begin(), syst.type.end(),
      [](const auto &entry) { return IsFusedOsc(entry.second); });
  if (!hasFloatingOscillation)
    for (const auto &dataset : setup.dsPDFMap)
      for (const auto &component : dataset.second)
        if (component.second.GetOscillated())
          throw std::runtime_error("MCMC needs a FusedOsc systematic for oscillated components");

  // Build osc pars
  OscParams osc = BuildOscParams(in.p.noms);

  // Prune any parameters that aren't in fit config
  DatasetParsInfo dsParsInfo = PruneUnusedParameters(in, setup, syst, osc);

  PrintFitSummary(in, dsParsInfo.datasetPars, osc);

  // Builds the components of the Asimov datsets with nominal systs applied, and combines with
  // data to make a llh for each dataset
  BuildResult datasets = BuildDatasets(in, dirs, setup, syst, osc, indexDistance, dsParsInfo);
  
  // Now combine the LLH for each dataset
  CombinedLLH combined = CombineLLHs(datasets.llhs, datasets.initialByDataset);
  StatisticSum &fullLLH = combined.llh;
  ParameterDict allParVals = combined.initialValues;

  // Finally apply all our constraints
  ApplyConstraints(fullLLH, in.constraints, in.p.noms);
  fullLLH.RegisterFitComponents();

  // Sampling settings remain specific to MCMC
  const int nsteps = fitConfig.GetIterations();
  const int burnin = fitConfig.GetBurnIn();
  const double sigmaScale = fitConfig.GetSigmaScale();
  const std::string &outDir = dirs.outDir;
  const std::string &postfitDistDir = dirs.postfitDistDir;
  ParameterDict &mins = in.p.mins;
  ParameterDict &maxs = in.p.maxs;
  ParameterDict &noms = in.p.noms;
  ParameterDict sigmas = in.p.sigmas;
  const ParameterDict &nbins = in.p.nbins;
  const ParameterDict &constrSigmas = in.constraints.sigmas;
  auto &fixedPars = in.p.fixedPars;

  // Starting point:
  //  - a fixed parameter stays nominal (the sampler never moves it);
  //  - the oscillation parameters (those of the FusedOsc systematics), and any free parameter with no Gaussian
  //    prior on itself, start uniformly in [min, max];
  //  - every other free parameter starts at nominal + Gaussian(0, prior sigma), redrawn until inside [min, max].
  std::set<std::string> oscParNames;
  for (std::map<std::string, std::string>::iterator systIt = syst.type.begin(); systIt != syst.type.end(); ++systIt)
    if (IsFusedOsc(systIt->second))
      oscParNames.insert(syst.paramNames[systIt->first].begin(), syst.paramNames[systIt->first].end());

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

  // And now the optimiser
  MetropolisSampler sampler;
  for (ParameterDict::iterator it = sigmas.begin(); it != sigmas.end(); ++it)
    sigmas[it->first] *= sigmaScale;
  sampler.SetSigmas(sigmas);
  for (std::map<std::string, bool>::iterator fixIt = fixedPars.begin(); fixIt != fixedPars.end(); ++fixIt)
    if (fixIt->second)
      sampler.Fix(fixIt->first);

  MCMC mh(sampler);
  mh.SetSaveChain(in.run.saveOutputs);
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

  // Create some axes for the mc to fill
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
  if (in.run.saveOutputs)
  {
    res.SaveAs(outDir + "/fit_result.txt");
    std::ofstream file(outDir + "/fit_result.txt", std::ios::app);
    file << "\nLLH: " << finalLLH << "\n";
    file << "\nFit Valid: " << validFit << std::endl;
    std::cout << "Saved fit result to " << outDir + "/fit_result.txt" << std::endl;
  }

  MCMCSamples samples = mh.GetSamples();

  if (in.run.saveOutputs)
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

  if (in.run.saveOutputs)
    std::cout << "Saving LH projections to \n\t" << projDir1D << "\n\t" << projDir2D << std::endl;
  std::cout << "Posterior (after burn-in, from the 1D projections): mean +- std" << std::endl;
  for (HistMap::const_iterator it = proj1D.begin(); it != proj1D.end(); ++it)
  {
    if (in.run.saveOutputs)
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
  if (in.run.saveOutputs)
    for (HistMap::const_iterator it = proj2D.begin(); it != proj2D.end(); ++it)
      IO::SaveHistogram(it->second, projDir2D + "/" + it->first + ".root");

  // Calculate postfit rates even when optional component spectra are not saved.
  RateByDatasetAndPdf bestFitNonBufferRates;
  for (DSMap::const_iterator dsIt = setup.dsPDFMap.begin(); dsIt != setup.dsPDFMap.end(); ++dsIt)
  {
    const DatasetModel &model = datasets.models.at(dsIt->first);
    BinnedED postfitDist = BinnedED("postfit dist", setup.systAxes);
    postfitDist.SetObservables(setup.dataObs);
    for (size_t i = 0; i < model.pdfs.size(); i++)
    {
      const std::string name = model.pdfs.at(i).GetName();
      BinnedED postfitPDF = ApplySystematicsByGroups(model.pdfs.at(i), model.pdfGroups.at(i),
                                                   syst.byDataset.at(dsIt->first), syst.group, bestFit);
      postfitPDF.Scale(bestFit.at(name));
      postfitPDF = ProjectToDataObs(postfitPDF, postfitDist, setup.dataObs);
      bestFitNonBufferRates[dsIt->first][name] += IntegrateBufferedBins(postfitPDF, setup.pdfConfig, setup.dataObs);
      if (in.run.saveOutputs)
        IO::SaveHistogram(postfitPDF.GetHistogram(), postfitDistDir + "/" + name + "_" + dsIt->first + ".root");
      // Sum all scaled distributions to get full postfit "dataset"
      postfitDist.Add(postfitPDF);
    }
    if (in.run.saveOutputs)
    {
      IO::SaveHistogram(postfitDist.GetHistogram(), postfitDistDir + "/postfitdist_" + dsIt->first + ".root");
      BinnedED dataDist = ProjectToDataObs(datasets.dataDists.at(dsIt->first), postfitDist, setup.dataObs);
      IO::SaveHistogram(dataDist.GetHistogram(), outDir + "/data_" + dsIt->first + ".root");
    }
  }
  PrintNominalBestFitSummary(in, bestFit, datasets.nominalNonBufferRates, bestFitNonBufferRates);

  if (in.run.saveOutputs)
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
