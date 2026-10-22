#include <FitParameters.hh>
#include <algorithm>
#include <functional>
#include <iostream>
#include <set>

namespace antinufit
{
  DatasetParameterMap PruneUnusedParameters(FitParameters &p,
                                            const DatasetSetup &setup,
                                            const SystConfig &systConfig,
                                            const OscParams &osc)
  {
    DatasetParameterMap datasetPars;

    ParameterDict &noms = p.noms;
    ParameterDict &mins = p.mins;
    ParameterDict &maxs = p.maxs;
    ParameterDict &sigmas = p.sigmas;
    std::map<std::string, bool> &fixedPars = p.fixedPars;
    const std::map<std::string, std::vector<std::string>> systParams = systConfig.GetParamNames();
    const std::map<std::string, std::vector<std::string>> systDatasets = systConfig.GetDataSets();

    // Iterate noms and prune
    for (ParameterDict::iterator it = noms.begin(); it != noms.end();)
    {
      const std::string &par = it->first;

      std::set<std::string> appliesTo;
      const bool isOsc = (par == "deltam21" || par == osc.theta12name);

      for (DSMap::const_iterator dsIt = setup.dsPDFMap.begin(); dsIt != setup.dsPDFMap.end(); ++dsIt)
      {
        const std::string &dsName = dsIt->first;
        const EvMap &evMap = dsIt->second;

        if (evMap.find(par) != evMap.end() || isOsc)
          appliesTo.insert(dsName);
      }

      for (const std::pair<const std::string, std::vector<std::string>> &entry : systParams)
        if (std::find(entry.second.begin(), entry.second.end(), par) != entry.second.end())
        {
          const std::vector<std::string> &datasets = systDatasets.at(entry.first);
          appliesTo.insert(datasets.begin(), datasets.end());
        }

      if (appliesTo.empty())
      {
        std::cout << par << " parameter defined in fit config but not in syst or event config. It will be ignored.\n";

        p.constraintSigmas.erase(par);
        p.constraintMeans.erase(par);

        p.ratioMeans.erase(par);
        p.ratioSigmas.erase(par);
        p.ratioParName.erase(par);

        p.fracMeans.erase(par);
        p.fracSigmas.erase(par);
        p.fracParName.erase(par);

        p.shapeMeans.erase(par);
        p.shapeSigmas.erase(par);
        p.shapeParNames.erase(par);
        p.shapeFuncName.erase(par);

        p.corrs.erase(par);
        p.corrParName.erase(par);

        mins.erase(par);
        maxs.erase(par);
        sigmas.erase(par);
        p.nbins.erase(par);
        fixedPars.erase(par);
        p.fdValues.erase(par);

        it = noms.erase(it);
        continue;
      }

      for (const std::string &dsName : appliesTo)
        datasetPars[dsName].push_back(par);

      ++it;
    }

    // Remove constraints that *reference* removed parameters
    // (e.g. ratioParName[p] == removed parameter, correlation partner removed, etc.)
    const std::function<bool(const std::string &)> isAlive = [&](const std::string &p)
    { return noms.find(p) != noms.end(); };

    // Ratio constraints: drop if either side missing
    for (std::map<std::string, std::string>::iterator it = p.ratioParName.begin(); it != p.ratioParName.end();)
    {
      const std::string &p1 = it->first;
      const std::string &p2 = it->second;
      if (!isAlive(p1) || !isAlive(p2))
      {
        p.ratioMeans.erase(p1);
        p.ratioSigmas.erase(p1);
        it = p.ratioParName.erase(it);
      }
      else
        ++it;
    }

    // Frac constraints: drop if either side missing
    for (std::map<std::string, std::string>::iterator it = p.fracParName.begin(); it != p.fracParName.end();)
    {
      const std::string &p1 = it->first;
      const std::string &p2 = it->second;
      if (!isAlive(p1) || !isAlive(p2))
      {
        p.fracMeans.erase(p1);
        p.fracSigmas.erase(p1);
        it = p.fracParName.erase(it);
      }
      else
        ++it;
    }

    // Correlations: drop if either side missing
    for (std::map<std::string, std::string>::iterator it = p.corrParName.begin(); it != p.corrParName.end();)
    {
      const std::string &p1 = it->first;
      const std::string &p2 = it->second;
      if (!isAlive(p1) || !isAlive(p2))
      {
        p.corrs.erase(p1);
        it = p.corrParName.erase(it);
      }
      else
        ++it;
    }

    return datasetPars;
  }
} // namespace antinufit
