#pragma once

#include <map>
#include <set>
#include <string>
#include <vector>

#include <SystSetup.hh>

#include <ParameterDict.h>
#include <BinnedED.h>

namespace antinufit
{

  // Apply all systematics that match any of `groupsForPdf` (or whose group == "").
  // Sets each systematic's parameters from `paramValues` before applying.
  // Returns the transformed distribution.
  inline BinnedED ApplySystematicsByGroups(
      BinnedED dist,
      const std::vector<std::string> &groupsForPdf,
      const std::map<std::string, Systematic *> &systsForDataset,
      const std::map<std::string, std::string> &systGroup,
      const ParameterDict &paramValues)
  {
    for (const std::string &grp : groupsForPdf)
    {
      for (std::map<std::string, Systematic *>::const_iterator it = systsForDataset.begin(); it != systsForDataset.end(); ++it)
      {
        const std::string &systName = it->first;
        Systematic *syst = it->second;

        // If group is "", apply to all groups
        const std::map<std::string, std::string>::const_iterator gIt = systGroup.find(systName);
        const std::string group = (gIt == systGroup.end()) ? "" : gIt->second;

        if (group == "" || group == grp)
        {
          // Set syst parameter values
          const std::set<std::string> pnames = syst->GetParameterNames();
          for (std::set<std::string>::const_iterator pit = pnames.begin(); pit != pnames.end(); ++pit)
          {
            const ParameterDict::const_iterator vIt = paramValues.find(*pit);
            if (vIt != paramValues.end())
              syst->SetParameter(*pit, vIt->second);
          }

          syst->Construct();
          double norm;
          dist = (*syst)(dist, &norm);
        }
      }
    }
    return dist;
  }

} // namespace antinufit
