#pragma once

#include <map>
#include <string>
#include <vector>

#include <BinnedED.h>
#include <BinnedNLLH.h>

namespace antinufit
{

  struct DatasetModel
  {
    // Model pieces for LH construction
    std::vector<BinnedED> pdfs;
    std::vector<std::vector<std::string>> pdfGroups;
    std::vector<int> genRates;
    std::vector<NormFittingStatus> normFittingStatuses;

    // Data products you want later
    BinnedED asimov;
    BinnedED fakeDataset;

    BinnedED dataDist;
  };

  using DatasetModelMap = std::map<std::string, DatasetModel>;

} // namespace antinufit
