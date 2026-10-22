#pragma once

#include <map>
#include <string>
#include <vector>

#include <BinnedNLLH.h>
#include <StatisticSum.h>
#include <TestStatistic.h>
#include <BuildDatasets.hh>

namespace antinufit
{

  // Use the first value found for parameters shared by multiple datasets
  inline ParameterDict GetInitialValues(const PreparedDatasetMap &datasets)
  {
    ParameterDict all;
    for (PreparedDatasetMap::const_iterator dataset = datasets.begin();
         dataset != datasets.end(); ++dataset)
    {
      for (ParameterDict::const_iterator p = dataset->second.initialValues.begin();
           p != dataset->second.initialValues.end(); ++p)
      {
        if (all.find(p->first) == all.end())
          all[p->first] = p->second;
      }
    }
    return all;
  }

  // Bring LLHs for each dataset together
  inline StatisticSum CombineLLHs(std::vector<BinnedNLLH> &testStats)
  {
    std::vector<TestStatistic *> ptrs;
    ptrs.reserve(testStats.size());
    for (std::vector<BinnedNLLH>::iterator lh = testStats.begin(); lh != testStats.end(); ++lh)
      ptrs.push_back(&(*lh));

    return Sum(ptrs);
  }

} // namespace antinufit
