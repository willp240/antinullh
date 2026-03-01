#pragma once

#include <map>
#include <string>
#include <vector>

#include <BinnedNLLH.h>
#include <StatisticSum.h>
#include <TestStatistic.h>

namespace antinufit
{

  struct CombinedLLH
  {
    StatisticSum llh;
    ParameterDict initialValues;
  };

  // Make parameterdict for all parameters for all datasets (without repeating)
  inline ParameterDict MergeInitialValues(const std::map<std::string, ParameterDict> &perDataset)
  {
    ParameterDict all;
    for (std::map<std::string, ParameterDict>::const_iterator kv = perDataset.begin();
         kv != perDataset.end(); ++kv)
    {
      for (ParameterDict::const_iterator p = kv->second.begin(); p != kv->second.end(); ++p)
      {
        if (all.find(p->first) == all.end())
          all[p->first] = p->second;
      }
    }
    return all;
  }

  // Bring LLHs for each dataset together
  inline CombinedLLH CombineLLHs(std::vector<BinnedNLLH> &testStats,
                                 const std::map<std::string, ParameterDict> &parameterValues)
  {
    std::vector<TestStatistic *> ptrs;
    ptrs.reserve(testStats.size());
    for (std::vector<BinnedNLLH>::iterator lh = testStats.begin(); lh != testStats.end(); ++lh)
      ptrs.push_back(&(*lh));

    CombinedLLH out{Sum(ptrs), MergeInitialValues(parameterValues)};
    return out;
  }

} // namespace antinufit
