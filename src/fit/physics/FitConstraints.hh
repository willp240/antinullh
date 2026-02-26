#pragma once

#include <map>
#include <string>
#include <vector>
#include <algorithm>

#include <StatisticSum.h>
#include <ShapeConstraints.hh>

namespace antinufit
{

  struct FitConstraints
  {
    // Simple Gaussian priors
    ParameterDict means;
    ParameterDict sigmas;

    // Correlations
    ParameterDict corrs;
    std::map<std::string, std::string> corrParName;

    // Ratio constraints
    ParameterDict ratioMeans;
    ParameterDict ratioSigmas;
    std::map<std::string, std::string> ratioParName;

    // Fractional diff constraints
    ParameterDict fracMeans;
    ParameterDict fracSigmas;
    std::map<std::string, std::string> fracParName;

    // Shape constraints
    ParameterDict shapeMeans;
    ParameterDict shapeSigmas;
    std::map<std::string, std::vector<std::string>> shapeParNames;
    std::map<std::string, std::string> shapeFuncName;
  };

  inline void ApplyConstraints(StatisticSum &fullLLH,
                               const FitConstraints &c,
                               const ParameterDict &noms)
  {
    // Correlations first
    std::vector<std::string> corrPairs;
    corrPairs.reserve(c.corrs.size());

    for (ParameterDict::const_iterator it = c.corrs.begin(); it != c.corrs.end(); ++it)
    {
      const std::string &p1 = it->first;
      const std::string &p2 = c.corrParName.at(p1);

      fullLLH.SetConstraint(p1,
                            c.means.at(p1), c.sigmas.at(p1),
                            p2,
                            c.means.at(p2), c.sigmas.at(p2),
                            it->second);

      corrPairs.push_back(p2);
    }

    // Single-parameter Gaussian priors (skip ones already covered by correlations)
    for (ParameterDict::const_iterator it = c.means.begin(); it != c.means.end(); ++it)
    {
      const std::string &p = it->first;
      if (c.corrs.find(p) != c.corrs.end())
        continue;
      if (std::find(corrPairs.begin(), corrPairs.end(), p) != corrPairs.end())
        continue;

      fullLLH.SetConstraint(p, it->second, c.sigmas.at(p));
    }

    // Ratio constraints
    for (ParameterDict::const_iterator it = c.ratioMeans.begin(); it != c.ratioMeans.end(); ++it)
      fullLLH.SetConstraint(it->first, c.ratioParName.at(it->first), it->second, c.ratioSigmas.at(it->first));

    // Fractional diff constraints
    for (ParameterDict::const_iterator it = c.fracMeans.begin(); it != c.fracMeans.end(); ++it)
      fullLLH.SetConstraint(it->first, c.fracParName.at(it->first), it->second, c.fracSigmas.at(it->first));

    // Shape constraints
    for (ParameterDict::const_iterator it = c.shapeMeans.begin(); it != c.shapeMeans.end(); ++it)
    {
      ShapeFunc func = GetShapeConstrFunc(c.shapeFuncName.at(it->first));
      fullLLH.SetConstraint(noms, func, it->second, c.shapeSigmas.at(it->first));
    }
  }

} // namespace antinufit
