// Fit the fixed-oscillation model with Minuit
// Reactor oscillation is baked into the PDFs at the configured parameter values
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
#include <RateReporting.hh>
#include <BuildLikelihood.hh>
#include <WriteDists.hh>

// OXO headers
#include <StatisticSum.h>
#include <IO.h>
#include <Rand.h>
#include <Minuit.h>

// ROOT headers
#include <TStopwatch.h>
#include <TMatrixD.h>
#include <TFile.h>

using namespace antinufit;

void fixedosc_fit(const std::string &fitConfigFile_,
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

  int numFixed = 0; // Will use this to count fixed pars later

  // Make sure output dirs exist
  OutputDirs dirs = MakeOutputDirs(fitConfig.GetOutDir());
  if (fitConfig.GetSaveOutputs())
    EnsureDirs(dirs);

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

  TStopwatch timer;
  timer.Start(true);
  fullLLH.Evaluate();
  timer.Stop();
  std::cout << "Eval time: " << timer.RealTime() << std::endl;

  // Now do a fit!
  Minuit min;
  min.SetMethod(fitConfig.GetMinuitMethod());
  min.SetMaxCalls(fitConfig.GetIterations());
  min.SetTolerance(fitConfig.GetMinuitTolerance());
  min.SetStrategy(fitConfig.GetMinuitStrategy());
  min.SetMinima(fitParameters.mins);
  min.SetMaxima(fitParameters.maxs);
  min.SetInitialValues(allParVals);
  min.SetInitialErrors(fitParameters.sigmas);
  for (std::map<std::string, bool>::iterator fixParMapIt = fitParameters.fixedPars.begin(); fixParMapIt != fitParameters.fixedPars.end(); fixParMapIt++)
  {
    if (fixParMapIt->second)
    {
      numFixed++;
      min.Fix(fixParMapIt->first);
    }
  }
  std::cout << "Run rabbit run!" << std::endl;
  TStopwatch minuitTimer;
  minuitTimer.Start(true);
  FitResult res = min.Optimise(&fullLLH);

  // MnMinimize may finish via its Simplex fallback.  That can locate a valid
  // minimum without producing a covariance matrix, which is not sufficient
  // for the detailed post-fit output.  In that mode, use the located point to
  // seed a final Migrad pass so Minuit calculates the Hessian/covariance
  if (fitConfig.GetMinuitMethod() == "Minimize")
  {
    const ParameterDict minimiseBestFit = res.GetBestFit();
    std::cout << "Minimize search complete; running Migrad from that point to obtain covariance." << std::endl;

    Minuit covarianceMin;
    covarianceMin.SetMethod("Migrad");
    covarianceMin.SetMaxCalls(fitConfig.GetIterations());
    covarianceMin.SetTolerance(fitConfig.GetMinuitTolerance());
    covarianceMin.SetStrategy(fitConfig.GetMinuitStrategy());
    covarianceMin.SetMinima(fitParameters.mins);
    covarianceMin.SetMaxima(fitParameters.maxs);
    covarianceMin.SetInitialValues(minimiseBestFit);
    covarianceMin.SetInitialErrors(fitParameters.sigmas);
    for (std::map<std::string, bool>::iterator fixParMapIt = fitParameters.fixedPars.begin(); fixParMapIt != fitParameters.fixedPars.end(); ++fixParMapIt)
    {
      if (fixParMapIt->second)
        covarianceMin.Fix(fixParMapIt->first);
    }
    res = covarianceMin.Optimise(&fullLLH);
  }

  minuitTimer.Stop();
  std::cout << "Minuit time: " << minuitTimer.RealTime() << std::endl;
  res.SetPrintPrecision(4);
  res.Print();
  ParameterDict bestFit = res.GetBestFit();
  fullLLH.SetParameters(bestFit);
  bool validFit = res.GetValid();
  double finalLLH = fullLLH.Evaluate();

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

  // Save the fit result and covariance matrix
  if (fitConfig.GetSaveOutputs())
  {
    res.SaveAs(dirs.outDir + "/fit_result.txt");
    std::ofstream file(dirs.outDir + "/fit_result.txt", std::ios::app);
    file << "\nLLH: " << finalLLH << "\n";
    file << "\nFit Valid: " << validFit << std::endl;
    file.close();
    TFile *outFile = new TFile((dirs.outDir + "/fit_result.root").c_str(), "RECREATE");
    DenseMatrix covMatrix = res.GetCovarianceMatrix();
    std::vector<std::string> paramNames;
    std::vector<std::string> allParamNames;
    std::vector<double> paramVals;
    std::vector<double> paramErr;
    TMatrixD covTMatrixD(bestFit.size() - numFixed, bestFit.size() - numFixed);
    for (ParameterDict::iterator parIt = bestFit.begin(); parIt != bestFit.end(); ++parIt)
    {
      allParamNames.push_back(parIt->first);

      // Fixed params won't be in the covariance matrix
      if (fitParameters.fixedPars[parIt->first])
      {
        continue;
      }

      paramNames.push_back(parIt->first);
      paramVals.push_back(bestFit[parIt->first]);
      if (validFit)
      {
        paramErr.push_back(sqrt(covMatrix.GetComponent(paramNames.size() - 1, paramNames.size() - 1)));
        for (int iParam = 0; iParam < paramNames.size(); iParam++)
        {
          covTMatrixD[paramNames.size() - 1][iParam] = covMatrix.GetComponent(paramNames.size() - 1, iParam);
          covTMatrixD[iParam][paramNames.size() - 1] = covMatrix.GetComponent(iParam, paramNames.size() - 1);
        }
      }
    }
    outFile->WriteObject(&paramNames, "paramNames");       // Non-fixed parameters
    outFile->WriteObject(&allParamNames, "allParamNames"); // All parameters
    outFile->WriteObject(&paramVals, "paramVals");
    outFile->WriteObject(&paramErr, "paramErr");
    outFile->WriteObject(&covTMatrixD, "covMatrix");
    std::cout << "Saved fit result to " << dirs.outDir + "/fit_result.txt and " << dirs.outDir << "/fit_result.root" << std::endl;
  }

  // Print some results to screen
  std::cout << "Fit complete for:" << std::endl;
  std::cout << "deltam: " << osc.deltam21 << std::endl;
  std::cout << osc.theta12name << ": " << osc.theta12 << std::endl;
  std::cout << "LLH: " << finalLLH << std::endl;
  for (PreparedDatasetMap::const_iterator dataset = datasets.begin(); dataset != datasets.end(); ++dataset)
    for (std::map<std::string, double>::const_iterator reacRatioIt = dataset->second.reactorRatio.begin();
         reacRatioIt != dataset->second.reactorRatio.end(); ++reacRatioIt)
      std::cout << reacRatioIt->first << " Reactor Ratio: " << reacRatioIt->second << std::endl;
  std::cout << "FitValid: " << validFit << std::endl;
  std::cout << std::endl
            << std::endl;
}

int main(int argc, char *argv[])
{
  if (argc != 5)
  {
    std::cout << "\nUsage fixedosc_fit <fit_config_file> <eve_config_file> <pdf_config_file> <syst_config_file>" << std::endl;
    return 1;
  }

  std::string fitConfigFile(argv[1]);
  std::string eveConfigFile(argv[2]);
  std::string pdfConfigPath(argv[3]);
  std::string systConfigFile(argv[4]);

  fixedosc_fit(fitConfigFile, eveConfigFile, pdfConfigPath, systConfigFile);

  return 0;
}
