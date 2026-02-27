// Antinu headers
#include <FitConfigLoader.hh>
#include <OscGridConfigLoader.hh>
#include <DatasetSetup.hh>
#include <OutputDirs.hh>
#include <FitSummary.hh>
#include <ApplySystematics.hh>
#include <Projection.hh>
#include <ParamUtils.hh>
#include <CombineLLHs.hh>
#include <BuildDatasets.hh>

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
                  const std::string &systConfigFile_,
                  const std::string &oscGridConfigFile_)
{
  Rand::SetSeed(0);

  // Load up the fit configuration information
  FitConfigLoader fitLoader(fitConfigFile_);
  FitConfig fitConfig = fitLoader.LoadActive();
  FitInputs in = MakeFitInputs(fitConfig);

  int numFixed = 0; // Will use this to count fixed pars later

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

  // Prune any parameters that aren't in fit config
  DatasetParsInfo dsParsInfo = PruneUnusedParameters(in, setup, syst, osc);

  PrintFitSummary(in, dsParsInfo.datasetPars, osc);

  // Builds the components of the Asimov datsets with nominal systs applied, and combines with
  // data to make a llh for each dataset
  BuildResult datasets = BuildDatasets(
      in, dirs, setup, syst, osc, indexDistance, dsParsInfo);

  // Now combine the LLH for each dataset
  CombinedLLH combined = CombineLLHs(datasets.llhs, datasets.initialByDataset);
  StatisticSum &fullLLH = combined.llh;
  ParameterDict allParVals = combined.initialValues;

  // Finally apply all our constraints
  ApplyConstraints(fullLLH, in.constraints, in.p.noms);
  fullLLH.RegisterFitComponents();

  TStopwatch timer;
  timer.Start(true);
  fullLLH.Evaluate();
  timer.Stop();
  std::cout << "Eval time: " << timer.RealTime() << std::endl;

  // Now do a fit!
  Minuit min;
  min.SetMethod(in.run.method);
  min.SetMaxCalls(in.run.iterations);
  min.SetTolerance(in.run.tolerance);
  min.SetStrategy(in.run.strategy);
  min.SetMinima(in.p.mins);
  min.SetMaxima(in.p.maxs);
  min.SetInitialValues(allParVals);
  min.SetInitialErrors(in.p.sigmas);
  for (std::map<std::string, bool>::iterator fixParMapIt = in.p.fixedPars.begin(); fixParMapIt != in.p.fixedPars.end(); fixParMapIt++)
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
  minuitTimer.Stop();
  std::cout << "Minuit time: " << minuitTimer.RealTime() << std::endl;
  res.SetPrintPrecision(4);
  res.Print();
  ParameterDict bestFit = res.GetBestFit();
  fullLLH.SetParameters(bestFit);
  bool validFit = res.GetValid();
  double finalLLH = fullLLH.Evaluate();

  // Now save the postfit distributions for each dataset
  for (DSMap::const_iterator dsIt = setup.dsPDFMap.begin(); dsIt != setup.dsPDFMap.end(); ++dsIt)
  {
    const std::string &dsName = dsIt->first;
    DatasetModel &m = datasets.models.at(dsName);

    BinnedED postfitDist("postfit dist", setup.systAxes);
    postfitDist.SetObservables(setup.dataObs);

    if (in.run.saveOutputs)
    {
      std::cout << "Saving scaled histograms and data for " << dsName
                << " to \n\t" << dirs.postfitDistDir << std::endl;
    }

    for (size_t i = 0; i < m.pdfs.size(); ++i)
    {
      const std::string pdfName = m.pdfs[i].GetName();

      BinnedED comp = m.pdfs[i];
      comp.Normalise();

      comp = ApplySystematicsByGroups(
          std::move(comp),
          m.pdfGroups[i],
          syst.byDataset.at(dsName),
          syst.group,
          bestFit);

      comp.Scale(bestFit.at(pdfName));

      // Project to whatever your "data obs" are (1D energy right now)
      BinnedED compToSave = ProjectToDataObs(comp, postfitDist, setup.dataObs);

      if (in.run.saveOutputs)
      {
        IO::SaveHistogram(compToSave.GetHistogram(),
                          dirs.postfitDistDir + "/" + pdfName + "_" + dsName + ".root");
      }

      postfitDist.Add(compToSave);
    }

    // Save the data with the same projection rule
    BinnedED dataToSave = ProjectToDataObs(datasets.dataDists.at(dsName), postfitDist, setup.dataObs);

    if (in.run.saveOutputs)
    {
      IO::SaveHistogram(postfitDist.GetHistogram(),
                        dirs.postfitDistDir + "/postfitdist_" + dsName + ".root");

      IO::SaveHistogram(dataToSave.GetHistogram(),
                        dirs.outDir + "/data_" + dsName + ".root");
    }
  }

  if (in.run.saveOutputs)
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
      if (in.p.fixedPars[parIt->first])
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

  std::cout << "Fit complete for:" << std::endl;
  std::cout << "deltam: " << osc.deltam21 << std::endl;
  std::cout << osc.theta12name << ": " << osc.theta12 << std::endl;
  std::cout << "LLH: " << finalLLH << std::endl;
  for (std::map<std::string, double>::const_iterator reacRatioIt = datasets.reactorRatio.begin();
       reacRatioIt != datasets.reactorRatio.end(); ++reacRatioIt)
  {
    std::cout << reacRatioIt->first << " Reactor Ratio: " << reacRatioIt->second << std::endl;
  }
  std::cout << "FitValid: " << validFit << std::endl;
  std::cout << std::endl
            << std::endl;
}

int main(int argc, char *argv[])
{
  if (argc != 6)
  {
    std::cout << "\nUsage fixedosc_fit <fit_config_file> <eve_config_file> <pdf_config_file> <syst_config_file> <oscgrid_config_file>" << std::endl;
    return 1;
  }

  std::string fitConfigFile(argv[1]);
  std::string eveConfigFile(argv[2]);
  std::string pdfConfigPath(argv[3]);
  std::string systConfigFile(argv[4]);
  std::string oscgridConfigFile(argv[5]);

  fixedosc_fit(fitConfigFile, eveConfigFile, pdfConfigPath, systConfigFile, oscgridConfigFile);

  return 0;
}
