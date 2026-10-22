#pragma once

#include <FitConfig.hh>
#include <SystConfig.hh>
#include <DatasetSetup.hh>

#include <unordered_map>

namespace antinufit
{
  // The four input files resolved into one description of an analysis
  // FitParameters is the mutable working copy used by a run
  struct AnalysisConfig
  {
    FitConfig fit;
    DatasetSetup datasets;
    SystConfig systematics;
    std::unordered_map<int, double> reactorBaselines;
  };

  AnalysisConfig LoadAnalysisConfig(const std::string &fitPath,
                                    const std::string &eventPath,
                                    const std::string &pdfPath,
                                    const std::string &systPath);
}
