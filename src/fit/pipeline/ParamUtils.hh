#pragma once

#include <string>
#include <map>
#include <vector>
#include <functional>

#include <ParameterDict.h>
#include <SystSetup.hh>

namespace antinufit
{
  // Initialise parameters to nominal for a dataset
  inline ParameterDict InitialiseDatasetParams(
      const ParameterDict &noms,
      const std::map<std::string, std::vector<std::string>> &datasetPars,
      const std::string &dsName)
  {
    ParameterDict values;
    const std::vector<std::string> &pars = datasetPars.at(dsName);

    for (std::vector<std::string>::const_reference p : pars)
      values[p] = noms.at(p);

    return values;
  }

  struct DatasetParsInfo
  {
    std::map<std::string, std::vector<std::string>> datasetPars;
  };

  // Mutates in.p.* and in.constraints.* to remove unused params, returns mapping param -> datasets it applies to
  inline DatasetParsInfo PruneUnusedParameters(FitInputs &in,
                                               const DatasetSetup &setup,
                                               const SystSetup &syst,
                                               const OscParams &osc)
  {
    DatasetParsInfo out;

    ParameterDict &noms = in.p.noms;
    ParameterDict &mins = in.p.mins;
    ParameterDict &maxs = in.p.maxs;
    ParameterDict &sigmas = in.p.sigmas;
    std::map<std::string, bool> &fixedPars = in.p.fixedPars;

    // Constraints bundle
    FitConstraints &c = in.constraints;

    // Iterate noms and prune
    for (ParameterDict::iterator it = noms.begin(); it != noms.end();)
    {
      const std::string &par = it->first;

      bool isPDF = false;
      bool isSyst = false;
      bool isOsc = false;

      for (DSMap::const_iterator dsIt = setup.dsPDFMap.begin(); dsIt != setup.dsPDFMap.end(); ++dsIt)
      {
        const std::string &dsName = dsIt->first;
        const EvMap &evMap = dsIt->second;

        if (evMap.find(par) != evMap.end())
        {
          isPDF = true;
          out.datasetPars[dsName].push_back(par);
        }
        else if (par == "deltam21" || par == osc.theta12name)
        {
          isOsc = true;
          out.datasetPars[dsName].push_back(par);
        }
      }

      if (!isPDF && !isOsc && syst.parDataSets.find(par) != syst.parDataSets.end())
      {
        isSyst = true;
        for (std::vector<std::string>::const_reference ds : syst.parDataSets.at(par))
          out.datasetPars[ds].push_back(par);
      }

      if (!isPDF && !isSyst && !isOsc)
      {
        std::cout << par << " parameter defined in fit config but not in syst or event config. It will be ignored.\n";

        c.sigmas.erase(par);
        c.means.erase(par);

        c.ratioMeans.erase(par);
        c.ratioSigmas.erase(par);
        c.ratioParName.erase(par);

        c.fracMeans.erase(par);
        c.fracSigmas.erase(par);
        c.fracParName.erase(par);

        c.shapeMeans.erase(par);
        c.shapeSigmas.erase(par);
        c.shapeParNames.erase(par);
        c.shapeFuncName.erase(par);

        c.corrs.erase(par);
        c.corrParName.erase(par);

        mins.erase(par);
        maxs.erase(par);
        sigmas.erase(par);
        fixedPars.erase(par);
        in.p.fdValues.erase(par);

        it = noms.erase(it);
        continue;
      }

      ++it;
    }

    // Remove constraints that *reference* removed parameters
    // (e.g. ratioParName[p] == removed parameter, correlation partner removed, etc.)
    const std::function<bool(const std::string &)> isAlive = [&](const std::string &p)
    { return noms.find(p) != noms.end(); };

    // Ratio constraints: drop if either side missing
    for (std::map<std::string, std::string>::iterator it = c.ratioParName.begin(); it != c.ratioParName.end();)
    {
      const std::string &p1 = it->first;
      const std::string &p2 = it->second;
      if (!isAlive(p1) || !isAlive(p2))
      {
        c.ratioMeans.erase(p1);
        c.ratioSigmas.erase(p1);
        it = c.ratioParName.erase(it);
      }
      else
        ++it;
    }

    // Frac constraints: drop if either side missing
    for (std::map<std::string, std::string>::iterator it = c.fracParName.begin(); it != c.fracParName.end();)
    {
      const std::string &p1 = it->first;
      const std::string &p2 = it->second;
      if (!isAlive(p1) || !isAlive(p2))
      {
        c.fracMeans.erase(p1);
        c.fracSigmas.erase(p1);
        it = c.fracParName.erase(it);
      }
      else
        ++it;
    }

    // Correlations: drop if either side missing
    for (std::map<std::string, std::string>::iterator it = c.corrParName.begin(); it != c.corrParName.end();)
    {
      const std::string &p1 = it->first;
      const std::string &p2 = it->second;
      if (!isAlive(p1) || !isAlive(p2))
      {
        c.corrs.erase(p1);
        it = c.corrParName.erase(it);
      }
      else
        ++it;
    }

    return out;
  }

} // namespace antinufit
