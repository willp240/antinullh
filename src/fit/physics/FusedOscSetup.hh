#pragma once

#include <SystSetup.hh>
#include <DatasetSetup.hh>
#include <OscWeight1D.hh>
#include <set>
#include <stdexcept>

namespace antinufit
{
  // A source belongs to one component in one dataset and must precede energy response.
  inline void ValidateFusedOscSetup(const DatasetSetup &setup, const SystSetup &syst)
  {
    for (const std::pair<const std::string, std::string> &entry : syst.group)
      if (entry.second.empty())
        throw std::runtime_error("Floating oscillation requires an explicit group for " + entry.first);

    for (const std::pair<const std::string, EvMap> &dataset : setup.dsPDFMap)
    {
      std::map<std::string, OscWeight1D *> byGroup;
      const std::map<std::string, Systematic *> &systematics = syst.byDataset.at(dataset.first);
      for (const std::pair<const std::string, Systematic *> &entry : systematics)
        if (OscWeight1D *osc = dynamic_cast<OscWeight1D *>(entry.second))
        {
          if (!byGroup.emplace(syst.group.at(entry.first), osc).second)
            throw std::runtime_error("Multiple FusedOsc systematics in the same group");
        }
      for (const std::pair<const std::string, Systematic *> &entry : systematics)
        if (byGroup.count(syst.group.at(entry.first)) &&
            !dynamic_cast<OscWeight1D *>(entry.second))
          throw std::runtime_error("FusedOsc must have its own group");

      std::set<OscWeight1D *> used;
      for (const std::pair<const std::string, EventConfig> &component : dataset.second)
      {
        const std::vector<std::string> groups = component.second.GetGroup();
        for (size_t i = 0; i < groups.size(); ++i)
          if (byGroup.count(groups[i]) && (!component.second.GetOscillated() || i != 0))
            throw std::runtime_error("FusedOsc must be first and apply only to an oscillated component");
        if (!component.second.GetOscillated())
          continue;
        if (groups.empty() || !byGroup.count(groups.front()) || component.second.GetNumDimensions() != 3)
          throw std::runtime_error("Oscillated component " + component.first + " needs a first FusedOsc group and three source dimensions");
        if (!used.insert(byGroup.at(groups.front())).second)
          throw std::runtime_error("Each FusedOsc source must belong to a separate component");
      }
      if (used.size() != byGroup.size())
        throw std::runtime_error("Unused FusedOsc systematic in dataset " + dataset.first);
    }
  }
}
