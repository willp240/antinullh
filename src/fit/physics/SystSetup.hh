#pragma once

#include <map>
#include <memory>
#include <string>
#include <vector>
#include <stdexcept>
#include <algorithm>
#include <cmath>
#include <unordered_map>

#include <SystConfigLoader.hh>
#include <SystFactory.hh>
#include <DistBuilder.hh>
#include <Systematic.h>
#include <DatasetSetup.hh>
#include <LoadDataset.hh>
#include <Event.h>

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
          throw std::runtime_error("Systematic " + systName + " uses parameter " + p +
                                   ", which is not defined in the fit config");
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

    for (const auto &dataset : dsPDFMap)
      out.byDataset[dataset.first];

    // Build the configured systematics for each dataset
    for (std::map<std::string, std::string>::const_iterator kv = out.type.begin();
         kv != out.type.end(); ++kv)
    {
      const std::string &systName = kv->first;

      for (const std::string &ds : out.dataSets.at(systName))
        if (!dsPDFMap.count(ds))
          throw std::runtime_error("Unknown dataset " + ds + " for systematic " + systName);

      out.fullParamNameVec.insert(out.fullParamNameVec.end(),
                                  out.paramNames[systName].begin(),
                                  out.paramNames[systName].end());

      // Each dataset owns its systematics; their parameters are still shared by name.
      for (const std::string &ds : out.dataSets.at(systName))
      {
        ShapeFunction shapeFunction;
        if (out.type.at(systName) == "Shape:ReactorOscProb")
        {
          const auto &params = out.paramNames.at(systName);
          const std::string &group = out.group.at(systName);
          if (params.size() != 2 || params.front() != "deltam21" || group.empty() ||
              out.distObs.at(systName) != std::vector<std::string>{"energy"} ||
              out.transObs.at(systName) != out.distObs.at(systName))
            throw std::runtime_error("Invalid reactor oscillation Shape configuration: " + systName);

          // This function needs the unoscillated MC of exactly one component.
          const EventConfig *component = nullptr;
          for (const auto &entry : dsPDFMap.at(ds))
          {
            const auto groups = entry.second.GetGroup();
            if (std::find(groups.begin(), groups.end(), group) == groups.end())
              continue;
            if (component || groups.front() != group || entry.second.GetNumDimensions() != 1 ||
                entry.second.GetOscillated() || entry.second.GetFlat())
              throw std::runtime_error("Reactor oscillation Shape needs one unoscillated 1D component with its group first: " + systName);
            component = &entry.second;
          }
          if (!component)
            throw std::runtime_error("No MC component for reactor oscillation Shape " + systName + " in " + ds);
          for (const auto &entry : out.group)
          {
            const auto &datasets = out.dataSets.at(entry.first);
            if (entry.first != systName && (entry.second.empty() || entry.second == group) &&
                std::find(datasets.begin(), datasets.end(), ds) != datasets.end())
              throw std::runtime_error("Reactor oscillation Shape needs its own group: " + systName);
          }

          auto mc = LoadMCDataset(component->GetPrunedPath());
          if (!mc)
            throw std::runtime_error("Cannot load MC for reactor oscillation Shape " + systName);
          const BinAxis reactorAxis = DistBuilder::BuildAxes(pdfConfig, 3).GetAxis(1);
          for (size_t i = 0; i < mc->GetNEntries(); ++i)
          {
            const double reactor = mc->GetEntry(i).GetDatum("reactorIndex");
            if (!std::isfinite(reactor) || reactor < reactorAxis.GetMin() || reactor >= reactorAxis.GetMax() ||
                std::fabs(reactor - std::round(reactor)) > 1e-6 || !indexDistance.count(int(std::lround(reactor))))
              throw std::runtime_error("Invalid reactor index or missing baseline in " + component->GetName() + " event " + std::to_string(i));
          }
          const BinnedED source = DistBuilder::Build(systName, 3, pdfConfig, mc.get());
          shapeFunction = MakeReactorOscillationShape(source, indexDistance, params[1]);
        }

        std::unique_ptr<Systematic> syst(SystFactory::New(
            systName, out.type.at(systName), out.paramNames.at(systName), noms, shapeFunction));
        syst->SetAxes(DistBuilder::BuildAxes(pdfConfig, out.distObs.at(systName).size()));
        syst->SetTransformationObs(out.transObs.at(systName));
        syst->SetDistributionObs(out.distObs.at(systName));
        syst->Construct();
        out.byDataset[ds][systName] = syst.get();
        out.owned.push_back(std::move(syst));
      }
    }

    // A parameter shared by several systematics applies to the union of their datasets.
    out.parDataSets.clear();
    for (const std::pair<const std::string, std::vector<std::string>> &entry : out.paramNames)
      for (const std::string &parameter : entry.second)
        for (const std::string &dataset : out.dataSets.at(entry.first))
        {
          std::vector<std::string> &datasets = out.parDataSets[parameter];
          if (std::find(datasets.begin(), datasets.end(), dataset) == datasets.end())
            datasets.push_back(dataset);
        }
    return out;
  }

} // namespace antinufit
