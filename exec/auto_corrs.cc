#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include <TBranch.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TH1D.h>
#include <TLeaf.h>
#include <TObjArray.h>
#include <TTree.h>
#include <TVirtualFFT.h>

namespace
{
Long64_t ReadCount(const std::string &text)
{
  std::size_t end = 0;
  const Long64_t value = std::stoll(text, &end);
  if (end != text.size() || value < 0)
    throw std::runtime_error("Expected a non-negative integer, got " + text);
  return value;
}

int Run(int argc, char *argv[])
{
  if (argc < 2)
    throw std::runtime_error("Usage: ./bin/auto_corrs chain.root [--burn-in STEPS] [--max-lag LAG] [--output FILE]");

  const std::string inputName = argv[1];
  std::string outputName = inputName;
  if (outputName.size() >= 5 && outputName.substr(outputName.size() - 5) == ".root")
    outputName.resize(outputName.size() - 5);
  outputName += "_autoCorr.root";
  Long64_t burnIn = 0;
  Long64_t maxLag = 1000;
  for (int i = 2; i < argc; ++i)
  {
    const std::string option = argv[i];
    if (i + 1 >= argc)
      throw std::runtime_error("Missing value for " + option);
    const std::string value = argv[++i];
    if (option == "--burn-in")
      burnIn = ReadCount(value);
    else if (option == "--max-lag")
      maxLag = ReadCount(value);
    else if (option == "--output")
      outputName = value;
    else
      throw std::runtime_error("Unknown option " + option);
  }
  if (outputName == inputName)
    throw std::runtime_error("Output must differ from input");

  std::unique_ptr<TFile> input(TFile::Open(inputName.c_str(), "READ"));
  if (!input || input->IsZombie())
    throw std::runtime_error("Cannot open input file " + inputName);
  TTree *tree = nullptr;
  input->GetObject("posteriors", tree);
  if (!tree || tree->GetEntries() == 0)
    throw std::runtime_error("Input must contain a non-empty posteriors tree");

  std::vector<std::string> names;
  TObjArray *branches = tree->GetListOfBranches();
  for (int i = 0; i < branches->GetEntries(); ++i)
  {
    TBranch *branch = static_cast<TBranch *>(branches->At(i));
    const std::string name = branch->GetName();
    if (name == "Step" || name == "StepTime" || name == "Accepted" || name == "LogL")
      continue;
    TLeaf *leaf = branch->GetLeaf(name.c_str());
    if (!leaf || leaf->GetLeafCount() || leaf->GetLenStatic() != 1 || std::string(leaf->GetTypeName()) != "Double_t")
      throw std::runtime_error("Expected a scalar double parameter branch: " + name);
    names.push_back(name);
  }
  if (names.empty())
    throw std::runtime_error("No parameter branches found");
  TLeaf *step = tree->GetLeaf("Step");
  if (!step)
    throw std::runtime_error("Missing Step branch");

  // Bind each branch once to stable storage, then copy the values for each entry
  std::vector<double> current(names.size());
  std::vector<std::vector<double> > samples(names.size());
  tree->SetBranchStatus("*", false);
  tree->SetBranchStatus("Step", true);
  for (std::size_t j = 0; j < names.size(); ++j)
  {
    tree->SetBranchStatus(names[j].c_str(), true);
    if (tree->SetBranchAddress(names[j].c_str(), &current[j]) < 0)
      throw std::runtime_error("Cannot bind parameter branch " + names[j]);
  }
  Long64_t previousStep = -1;
  for (Long64_t i = 0; i < tree->GetEntries(); ++i)
  {
    if (tree->GetEntry(i) <= 0)
      throw std::runtime_error("Cannot read entry " + std::to_string(i));
    const Long64_t stepNumber = step->GetValueLong64();
    if (stepNumber <= previousStep)
      throw std::runtime_error("Step resets within the tree: run autocorrelations separately for each chain, before hadd");
    previousStep = stepNumber;
    if (stepNumber < burnIn)
      continue;
    for (std::size_t j = 0; j < names.size(); ++j)
    {
      if (!std::isfinite(current[j]))
        throw std::runtime_error("Non-finite sample for " + names[j]);
      samples[j].push_back(current[j]);
    }
  }
  tree->ResetBranchAddresses();
  input.reset();
  const std::size_t count = samples[0].size();
  if (count < 2)
    throw std::runtime_error("Need at least two samples after the burn-in cut");
  maxLag = std::min(maxLag, static_cast<Long64_t>(count - 1));
  std::cout << "Read " << count << " samples after Step >= " << burnIn << std::endl;

  // Zero padding prevents the FFT from correlating the end of the chain with its start
  int fftSize = 1;
  while (static_cast<std::size_t>(fftSize) < 2 * count)
  {
    if (fftSize > std::numeric_limits<int>::max() / 2)
      throw std::runtime_error("Chain too long for ROOT FFT");
    fftSize *= 2;
  }
  std::unique_ptr<TVirtualFFT> forward(TVirtualFFT::FFT(1, &fftSize, "R2C ES K"));
  std::unique_ptr<TVirtualFFT> inverse(TVirtualFFT::FFT(1, &fftSize, "C2R ES K"));
  if (!forward || !inverse)
    throw std::runtime_error("ROOT FFT support is unavailable");
  std::vector<double> centered(fftSize, 0.0);
  std::vector<double> real(fftSize), imaginary(fftSize), covariance(fftSize);
  std::unique_ptr<TFile> output(TFile::Open(outputName.c_str(), "RECREATE"));
  if (!output || output->IsZombie())
    throw std::runtime_error("Cannot create output file " + outputName);
  TDirectory *directory = output->mkdir("Auto_Corrs");
  if (!directory)
    throw std::runtime_error("Cannot create Auto_Corrs directory");
  directory->cd();
  for (std::size_t j = 0; j < names.size(); ++j)
  {
    const std::vector<double> &values = samples[j];
    if (*std::min_element(values.begin(), values.end()) == *std::max_element(values.begin(), values.end()))
    {
      std::cout << "Skipping constant parameter " << names[j] << std::endl;
      continue;
    }
    // Subtract a reference before averaging to retain precision for nearly constant values
    double meanOffset = 0.0;
    for (std::size_t i = 0; i < count; ++i)
      meanOffset += values[i] - values[0];
    meanOffset /= count;
    std::fill(centered.begin(), centered.end(), 0.0);
    for (std::size_t i = 0; i < count; ++i)
      centered[i] = (values[i] - values[0]) - meanOffset;
    forward->SetPoints(centered.data());
    forward->Transform();
    forward->GetPointsComplex(real.data(), imaginary.data());
    for (int k = 0; k <= fftSize / 2; ++k)
    {
      real[k] = real[k] * real[k] + imaginary[k] * imaginary[k];
      imaginary[k] = 0.0;
    }
    inverse->SetPointsComplex(real.data(), imaginary.data());
    inverse->Transform();
    inverse->GetPoints(covariance.data());
    if (!std::isfinite(covariance[0]) || covariance[0] <= 0.0)
      throw std::runtime_error("Invalid variance for " + names[j]);
    TH1D histogram(names[j].c_str(), names[j].c_str(), static_cast<int>(maxLag + 1), -0.5, maxLag + 0.5);
    histogram.GetXaxis()->SetTitle("Lag (steps)");
    histogram.GetYaxis()->SetTitle("Auto-correlation");
    for (Long64_t k = 0; k <= maxLag; ++k)
      histogram.SetBinContent(static_cast<int>(k + 1), covariance[k] / covariance[0]);
    if (histogram.Write() <= 0)
      throw std::runtime_error("Cannot write autocorrelation for " + names[j]);
  }
  output->Close();
  if (output->TestBit(TFile::kWriteError))
    throw std::runtime_error("Failed to write output file " + outputName);
  std::cout << "Saved " << outputName << " with lags 0 to " << maxLag << std::endl;
  return 0;
}
}

int main(int argc, char *argv[])
{
  try
  {
    return Run(argc, argv);
  }
  catch (const std::exception &error)
  {
    std::cerr << "auto_corrs: " << error.what() << std::endl;
    return 1;
  }
}
