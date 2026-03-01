#pragma once

#include <map>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

#include <SystConfigLoader.hh>
#include <SystFactory.hh>
#include <DistBuilder.hh>
#include <Systematic.h>

namespace antinufit
{

  struct SystSetup
  {
    // Own all Systematic objects
    std::vector<std::unique_ptr<Systematic>> owned;

    // Dataset -> (systName -> pointer)
    std::map<std::string, std::map<std::string, Systematic *>> byDataset;

    std::map<std::string, std::vector<std::string>> paramNames;  // systName -> [param...]
    std::map<std::string, std::string> group;                    // systName -> group
    std::map<std::string, std::string> type;                     // systName -> type
    std::map<std::string, std::vector<std::string>> distObs;     // systName -> dist obs
    std::map<std::string, std::vector<std::string>> transObs;    // systName -> trans obs
    std::map<std::string, std::vector<std::string>> dataSets;    // systName -> datasets
    std::map<std::string, std::vector<std::string>> parDataSets; // paramName -> datasets

    std::vector<std::string> fullParamNameVec;
  };

  inline void CheckSystParams(const SystSetup &s,
                              const ParameterDict &noms)
  {
    // Check any parameter defined for a systematic was declared in the fit config
    for (std::map<std::string, std::vector<std::string>>::const_iterator kv = s.paramNames.begin();
         kv != s.paramNames.end(); ++kv)
    {
      const std::string &systName = kv->first;
      const std::vector<std::string> &params = kv->second;
      (void)systName;

      for (std::vector<std::string>::const_reference p : params)
      {
        if (noms.find(p) == noms.end())
        {
          std::cout << "ERROR: Syst config defines systematic with parameter " << p
                    << ", but this parameter is not defined in the fit config" << std::endl;
          throw;
        }
      }
    }
  }

  inline SystSetup BuildSystematics(const std::string &systConfigFile,
                                    const PDFConfig &pdfConfig,
                                    const ParameterDict &noms,
                                    const DSMap &dsPDFMap,
                                    const std::unordered_map<int, double> &indexDistance)
  {
    SystSetup out;

    // Load config
    SystConfigLoader loader(systConfigFile);
    SystConfig cfg = loader.LoadActive();

    out.paramNames = cfg.GetParamNames();
    out.group = cfg.GetGroup();
    out.type = cfg.GetType();
    out.distObs = cfg.GetDistObs();
    out.transObs = cfg.GetTransObs();
    out.dataSets = cfg.GetDataSets();
    out.parDataSets = cfg.GetParDataSets();

    // Validate params exist in fit config
    CheckSystParams(out, noms);

    // Dummy osc grid map
    std::map<int, OscGrid *> dummyOscGridMap;

    // Build each systematic once
    for (std::map<std::string, std::string>::const_iterator kv = out.type.begin();
         kv != out.type.end(); ++kv)
    {
      const std::string &systName = kv->first;

      out.fullParamNameVec.insert(out.fullParamNameVec.end(),
                                  out.paramNames[systName].begin(),
                                  out.paramNames[systName].end());

      std::unique_ptr<Systematic> syst(
          SystFactory::New(systName,
                           out.type[systName],
                           out.paramNames[systName],
                           noms,
                           dummyOscGridMap,
                           indexDistance));

      AxisCollection axes = DistBuilder::BuildAxes(pdfConfig, out.distObs[systName].size());
      syst->SetAxes(axes);
      syst->SetTransformationObs(out.transObs[systName]);
      syst->SetDistributionObs(out.distObs[systName]);
      syst->Construct();

      // Assign to datasets
      for (std::vector<std::string>::const_reference ds : out.dataSets[systName])
      {
        if (dsPDFMap.find(ds) == dsPDFMap.end())
        {
          std::cout << "ERROR: SystConfig defines systematic " << systName
                    << " applying to dataset " << ds
                    << " that is not defined in EventConfig" << std::endl;
          throw;
        }
        out.byDataset[ds][systName] = syst.get();
      }

      out.owned.push_back(std::move(syst));
    }

    return out;
  }

} // namespace antinufit
