#include <BuildSysts.hh>
#include <SystConfigLoader.hh>
#include <SystFactory.hh>
#include <BuildDist.hh>
#include <LoadEvents.hh>
#include <SystFuncs.hh>
#include <Event.h>

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace antinufit
{
  namespace
  {
  void CheckSystParams(const std::map<std::string, std::vector<std::string>> &paramNames,
                       const ParameterDict &noms)
  {
    // Check any parameter defined for a systematic was declared in the fit config
    for (std::map<std::string, std::vector<std::string>>::const_iterator kv = paramNames.begin();
         kv != paramNames.end(); ++kv)
    {
      const std::string &systName = kv->first;
      const std::vector<std::string> &params = kv->second;
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

  } // namespace

  SystsByDataset BuildSysts(const SystConfig &cfg,
                                      const PDFConfig &pdfConfig,
                                      const ParameterDict &noms,
                                      const DSMap &dsPDFMap,
                                      const std::unordered_map<int, double> &indexDistance)
  {
    SystsByDataset out;

    const std::map<std::string, std::vector<std::string>> paramNames = cfg.GetParamNames();
    const std::map<std::string, std::string> groups = cfg.GetGroup();
    const std::map<std::string, std::string> type = cfg.GetType();
    const std::map<std::string, std::vector<std::string>> distObs = cfg.GetDistObs();
    const std::map<std::string, std::vector<std::string>> transObs = cfg.GetTransObs();
    const std::map<std::string, std::vector<std::string>> dataSets = cfg.GetDataSets();

    // Validate params exist in fit config
    CheckSystParams(paramNames, noms);

    for (const std::pair<const std::string, EvMap> &dataset : dsPDFMap)
      out[dataset.first];

    // Build the configured systematics for each dataset
    for (std::map<std::string, std::string>::const_iterator kv = type.begin();
         kv != type.end(); ++kv)
    {
      const std::string &systName = kv->first;

      for (const std::string &ds : dataSets.at(systName))
        if (!dsPDFMap.count(ds))
          throw std::runtime_error("Unknown dataset " + ds + " for systematic " + systName);

      // Create a separate object for each dataset while sharing parameter names
      for (const std::string &ds : dataSets.at(systName))
      {
        ShapeFunction shapeFunction;
        // This Shape function needs a detailed MC histogram before it can be created
        if (type.at(systName) == "Shape:ReactorOscProb")
        {
          // Check the parameter names, energy observables and group are configured correctly
          const std::vector<std::string> &params = paramNames.at(systName);
          const std::string &group = groups.at(systName);
          if (params.size() != 2 || params.front() != "deltam21" || group.empty() ||
              distObs.at(systName) != std::vector<std::string>{"energy"} ||
              transObs.at(systName) != distObs.at(systName))
            throw std::runtime_error("Invalid reactor oscillation Shape configuration: " + systName);

          // Find one matching component with an unoscillated 1D PDF and this group first
          const EventConfig *component = nullptr;
          for (const std::pair<const std::string, EventConfig> &entry : dsPDFMap.at(ds))
          {
            const std::vector<std::string> groups = entry.second.GetGroup();
            if (std::find(groups.begin(), groups.end(), group) == groups.end())
              continue;
            if (component || groups.front() != group || entry.second.GetNumDimensions() != 1 ||
                entry.second.GetOscillated() || entry.second.GetFlat())
              throw std::runtime_error("Reactor oscillation Shape needs one unoscillated 1D component with its group first: " + systName);
            component = &entry.second;
          }
          // Stop if no component supplies the MC for this Shape
          if (!component)
            throw std::runtime_error("No MC component for reactor oscillation Shape " + systName + " in " + ds);

          // Keep this group exclusive so another systematic cannot run before oscillation
          for (const std::pair<const std::string, std::string> &entry : groups)
          {
            const std::vector<std::string> &datasets = dataSets.at(entry.first);
            if (entry.first != systName && (entry.second.empty() || entry.second == group) &&
                std::find(datasets.begin(), datasets.end(), ds) != datasets.end())
              throw std::runtime_error("Reactor oscillation Shape needs its own group: " + systName);
          }

          // Open the selected component's MC file
          std::unique_ptr<DataSet> mc = LoadEvents(component->GetPrunedPath());
          if (!mc)
            throw std::runtime_error("Cannot load MC for reactor oscillation Shape " + systName);

          // Check reactor IDs and baselines before filling can clamp invalid IDs into edge bins
          const BinAxis reactorAxis = BuildDist::BuildAxes(pdfConfig, 3).GetAxis(1);
          for (size_t i = 0; i < mc->GetNEntries(); ++i)
          {
            const double reactor = mc->GetEntry(i).GetDatum("reactorIndex");
            if (!std::isfinite(reactor) || reactor < reactorAxis.GetMin() || reactor >= reactorAxis.GetMax() ||
                std::fabs(reactor - std::round(reactor)) > 1e-6 || !indexDistance.count(int(std::lround(reactor))))
              throw std::runtime_error("Invalid reactor index or missing baseline in " + component->GetName() + " event " + std::to_string(i));
          }

          // Build (reconstructed energy, reactor ID, true energy) counts and prepare the Shape function
          const BinnedED source = BuildDist::Build(systName, 3, pdfConfig, mc.get());
          shapeFunction = MakeReactorOscillationShape(source, indexDistance, params[1]);
        }

        std::unique_ptr<Systematic> syst(SystFactory::New(
            systName, type.at(systName), paramNames.at(systName), noms, shapeFunction));
        syst->SetAxes(BuildDist::BuildAxes(pdfConfig, distObs.at(systName).size()));
        syst->SetTransformationObs(transObs.at(systName));
        syst->SetDistributionObs(distObs.at(systName));
        syst->Construct();
        out[ds][systName] = std::move(syst);
      }
    }

    return out;
  }

} // namespace antinufit
