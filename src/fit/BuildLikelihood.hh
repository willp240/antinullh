#pragma once

#include <BuildDatasets.hh>

namespace antinufit
{
  void ConfigureDatasetLikelihood(
      BinnedNLLH &lh,
      const SystsByDataset &syst,
      const SystConfig &systConfig,
      const FitConfig &fitConfig,
      const PDFConfig &pdfConfig,
      const std::string &dsName,
      const PreparedDataset &dataset,
      const std::vector<BinnedED> &pdfs,
      const std::vector<int> &generatedRates);

  // Construct the OXO likelihoods from already prepared data and component PDFs
  // The returned vector must outlive any StatisticSum built from it
  std::vector<BinnedNLLH> BuildDatasetLikelihoods(
      const PreparedDatasetMap &datasets,
      const SystsByDataset &syst,
      const SystConfig &systConfig,
      const FitConfig &fitConfig,
      const PDFConfig &pdfConfig);

}
