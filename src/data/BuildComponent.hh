#pragma once

#include <BuildSysts.hh>

#include <BinnedED.h>
#include <ParameterDict.h>

#include <string>
#include <vector>

namespace antinufit
{
  // Distort and scale a component distribution, then project it onto the data observables
  BinnedED BuildComponent(
      const BinnedED &pdf,
      const std::vector<std::string> &pdfGroups,
      const SystsByDataset &syst,
      const std::map<std::string, std::string> &systGroups,
      const std::string &datasetName,
      const ParameterDict &values,
      double rateScale,
      const BinnedED &dataAxes,
      const std::vector<std::string> &dataObs);

  // Build each component at the supplied values and scales, adding it to total
  std::vector<BinnedED> BuildComponents(
      const std::vector<BinnedED> &pdfs,
      const std::vector<std::vector<std::string>> &pdfGroups,
      const SystsByDataset &syst,
      const SystConfig &systConfig,
      const std::string &datasetName,
      const ParameterDict &values,
      const std::vector<double> &rateScales,
      BinnedED &total,
      const std::vector<std::string> &dataObs);
}
