#pragma once

#include <map>
#include <string>
#include <vector>
#include <algorithm>

#include <StatisticSum.h>
#include <ShapeConstraints.hh>
#include <FitParameters.hh>

namespace antinufit
{

  inline void ApplyConstraints(StatisticSum &fullLLH,
                               const FitParameters &p,
                               const ParameterDict &noms)
  {
    // Correlations first
    std::vector<std::string> corrPairs;
    corrPairs.reserve(p.corrs.size());

    for (ParameterDict::const_iterator it = p.corrs.begin(); it != p.corrs.end(); ++it)
    {
      const std::string &p1 = it->first;
      const std::string &p2 = p.corrParName.at(p1);

      fullLLH.SetConstraint(p1,
                            p.constraintMeans.at(p1), p.constraintSigmas.at(p1),
                            p2,
                            p.constraintMeans.at(p2), p.constraintSigmas.at(p2),
                            it->second);

      corrPairs.push_back(p2);
    }

    // Single-parameter Gaussian priors (skip ones already covered by correlations)
    for (ParameterDict::const_iterator it = p.constraintMeans.begin(); it != p.constraintMeans.end(); ++it)
    {
      const std::string &parName = it->first;
      if (p.corrs.find(parName) != p.corrs.end())
        continue;
      if (std::find(corrPairs.begin(), corrPairs.end(), parName) != corrPairs.end())
        continue;

      fullLLH.SetConstraint(parName, it->second, p.constraintSigmas.at(parName));
    }

    // Ratio constraints
    for (ParameterDict::const_iterator it = p.ratioMeans.begin(); it != p.ratioMeans.end(); ++it)
      fullLLH.SetConstraint(it->first, p.ratioParName.at(it->first), it->second, p.ratioSigmas.at(it->first));

    // Fractional diff constraints
    for (ParameterDict::const_iterator it = p.fracMeans.begin(); it != p.fracMeans.end(); ++it)
      fullLLH.SetConstraint(it->first, p.fracParName.at(it->first), it->second, p.fracSigmas.at(it->first));

    // Shape constraints
    for (ParameterDict::const_iterator it = p.shapeMeans.begin(); it != p.shapeMeans.end(); ++it)
    {
      ShapeFunc func = GetShapeConstrFunc(p.shapeFuncName.at(it->first));
      fullLLH.SetConstraint(noms, func, it->second, p.shapeSigmas.at(it->first));
    }
  }

} // namespace antinufit
