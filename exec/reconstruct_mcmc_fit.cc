// Reconstruct non-buffer component yields and the selected best-sample spectra
#include <AnalysisConfig.hh>
#include <BuildComponent.hh>
#include <BuildDatasets.hh>
#include <FitParameters.hh>
#include <OscParams.hh>
#include <RateReporting.hh>
#include <WriteDists.hh>
#include <TFile.h>
#include <TTree.h>
#include <TLeaf.h>
#include <TObjArray.h>
#include <TStopwatch.h>
#include <nlohmann/json.hpp>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <memory>
#include <stdexcept>

using namespace antinufit;

int Run(int argc, char **argv)
{
  if (argc != 7)
    throw std::runtime_error("Usage: reconstruct_mcmc_fit fit.ini event.ini pdf.ini syst.ini requests.root output_directory");
  const std::string outputPath = std::filesystem::absolute(argv[6]).string();
  std::filesystem::create_directories(outputPath);
  TStopwatch setupTimer;
  setupTimer.Start();
  AnalysisConfig config = LoadAnalysisConfig(argv[1], argv[2], argv[3], argv[4]);
  config.fit.SetSaveOutputs(true);
  FitParameters parameters = MakeFitParameters(config.fit);
  const DatasetSetup &setup = config.datasets;
  SystsByDataset systematics = BuildSysts(config.systematics, setup.pdfConfig,
      parameters.noms, setup.dsPDFMap, config.reactorBaselines);
  OscParams osc = BuildOscParams(parameters.noms);
  DatasetParameterMap used = PruneUnusedParameters(parameters, setup, config.systematics, osc);
  PreparedDatasetMap datasets = BuildDatasets(config.fit, parameters, setup, systematics,
      config.systematics, osc, config.reactorBaselines, used, config.randomSeed);
  ApplyReactorRateScales(parameters, datasets);
  setupTimer.Stop();

  std::unique_ptr<TFile> input(TFile::Open(argv[5], "READ"));
  TTree *requests = nullptr;
  if (!input || input->IsZombie())
    throw std::runtime_error("Cannot open rate requests");
  input->GetObject("requests", requests);
  if (!requests || requests->GetEntries() == 0)
    throw std::runtime_error("Missing non-empty requests tree");
  ParameterDict values = parameters.noms;
  Long64_t sampleID = 0;
  if (requests->SetBranchAddress("sample_id", &sampleID) < 0)
    throw std::runtime_error("Missing sample_id branch");
  TObjArray *branches = requests->GetListOfBranches();
  for (int i = 0; i < branches->GetEntries(); ++i)
  {
    const std::string name = branches->At(i)->GetName();
    if (name == "sample_id")
      continue;
    if (!values.count(name))
      throw std::runtime_error("Unknown requested parameter " + name);
    if (requests->SetBranchAddress(name.c_str(), &values.at(name)) < 0)
      throw std::runtime_error("Invalid requested parameter " + name);
  }

  const OutputDirs dirs = MakeOutputDirs(outputPath + "/bestfit");
  EnsureDirs(dirs);
  WriteInputDists(datasets, setup, dirs);
  std::unique_ptr<TFile> output(TFile::Open((outputPath + "/rates.root").c_str(), "RECREATE"));
  if (!output || output->IsZombie())
    throw std::runtime_error("Cannot create rates.root");
  TTree rates("rates", "Expected component counts outside likelihood buffer bins");
  rates.Branch("sample_id", &sampleID, "sample_id/L");
  std::map<std::string, double> rateValues;
  nlohmann::json columns;
  for (PreparedDatasetMap::const_iterator dataset = datasets.begin(); dataset != datasets.end(); ++dataset)
    for (const BinnedED &pdf : dataset->second.pdfs)
    {
      const std::string name = "yield__" + dataset->first + "__" + pdf.GetName();
      rateValues[name] = 0.0;
      columns[name] = {{"dataset", dataset->first}, {"component", pdf.GetName()},
                       {"nominal", dataset->second.nominalNonBufferRates.at(pdf.GetName())},
                       {"units", "expected events outside configured buffer bins"}};
    }
  for (std::map<std::string, double>::iterator entry = rateValues.begin(); entry != rateValues.end(); ++entry)
    rates.Branch(entry->first.c_str(), &entry->second, (entry->first + "/D").c_str());

  TStopwatch evaluationTimer;
  evaluationTimer.Start();
  bool wroteBest = false;
  for (Long64_t i = 0; i < requests->GetEntries(); ++i)
  {
    if (requests->GetEntry(i) <= 0)
      throw std::runtime_error("Cannot read requested sample");
    for (PreparedDatasetMap::const_iterator dataset = datasets.begin(); dataset != datasets.end(); ++dataset)
    {
      BinnedED total("postfit dist", setup.systAxes);
      total.SetObservables(setup.dataObs);
      const std::vector<BinnedED> components = BuildComponents(dataset->second.pdfs,
          dataset->second.pdfGroups, systematics, config.systematics, dataset->first,
          values, std::vector<double>(dataset->second.pdfs.size(), 1.0), total, setup.dataObs);
      const std::map<std::string, double> yields = GetComponentRates(dataset->second.pdfs,
          components, setup.pdfConfig, setup.dataObs);
      for (std::map<std::string, double>::const_iterator entry = yields.begin(); entry != yields.end(); ++entry)
        rateValues.at("yield__" + dataset->first + "__" + entry->first) = entry->second;
      if (sampleID == -1)
        WritePostfitDists(dataset->first, dataset->second, setup, components, total, dirs);
    }
    if (sampleID == -1)
      wroteBest = true;
    rates.Fill();
  }
  evaluationTimer.Stop();
  requests->ResetBranchAddresses();
  if (!wroteBest)
    throw std::runtime_error("Requests must include sample_id = -1 for the selected best sample");
  output->cd();
  if (rates.Write() <= 0)
    throw std::runtime_error("Cannot write derived yields");
  // Detach the stack-owned tree before closing the ROOT file
  rates.SetDirectory(nullptr);
  output->Close();
  if (output->TestBit(TFile::kWriteError))
    throw std::runtime_error("Error writing derived yields");
  nlohmann::json metadata = {{"columns", columns}, {"evaluations", requests->GetEntries()},
      {"setup_seconds", setupTimer.RealTime()}, {"evaluation_seconds", evaluationTimer.RealTime()},
      {"seconds_per_evaluation", evaluationTimer.RealTime() / requests->GetEntries()},
      {"method", "BuildComponents and GetComponentRates at each requested posterior sample"}};
  std::ofstream description(outputPath + "/rate_metadata.json");
  description << metadata.dump(2) << std::endl;
  std::cout << metadata.dump(2) << std::endl;
  return 0;
}

int main(int argc, char **argv)
{
  try { return Run(argc, argv); }
  catch (const std::exception &error)
  {
    std::cerr << "reconstruct_mcmc_fit: " << error.what() << std::endl;
    return 1;
  }
}
