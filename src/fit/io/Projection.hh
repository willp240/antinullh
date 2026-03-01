#pragma once

#include <vector>
#include <BinnedED.h>

namespace antinufit
{

  // Ensure component matches the dimensionality of the data distribution.
  inline BinnedED ProjectToDataObs(const BinnedED &src,
                                   const BinnedED &dataDist,
                                   const std::vector<std::string> &dataObs)
  {
    if (src.GetNDims() != dataDist.GetNDims())
      return src.Marginalise(dataObs);

    return src;
  }

} // namespace antinufit
