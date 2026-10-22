#include <AnalysisConfig.hh>

#include <FitConfigLoader.hh>
#include <SystConfigLoader.hh>
#include <ReactorBaselines.hh>

#include <algorithm>
#include <set>
#include <stdexcept>

namespace antinufit
{
  AnalysisConfig LoadAnalysisConfig(const std::string &fitPath,
                                    const std::string &eventPath,
                                    const std::string &pdfPath,
                                    const std::string &systPath)
  {
    FitConfigLoader fitLoader(fitPath);
    SystConfigLoader systLoader(systPath);

    AnalysisConfig result;
    result.fit = fitLoader.LoadActive();
    result.datasets = LoadDatasetSetup(eventPath, pdfPath);
    result.systematics = systLoader.LoadActive();
    result.reactorBaselines = LoadIndexDistanceMap(result.fit.GetReactorsJsonFile());

    const ParameterDict parameters = result.fit.GetNominals();
    const DSMap &datasets = result.datasets.dsPDFMap;
    if (datasets.empty())
      throw std::runtime_error("No active datasets in event config");

    std::set<std::string> fixedReactorRates;
    for (DSMap::const_iterator dataset = datasets.begin(); dataset != datasets.end(); ++dataset)
    {
      if (dataset->second.empty())
        throw std::runtime_error("No active components in dataset " + dataset->first);
      for (EvMap::const_iterator component = dataset->second.begin();
           component != dataset->second.end(); ++component)
      {
        if (!parameters.count(component->first))
          throw std::runtime_error("Component " + component->first + " in " + dataset->first +
                                   " has no rate parameter in fit config");
        if (component->second.GetOscillated() &&
            !fixedReactorRates.insert(component->first).second)
          throw std::runtime_error("Fixed reactor component " + component->first +
                                   " is used in more than one dataset; its rate scaling is ambiguous");
      }
    }

    // Catch dataset names and parameter names that cannot be resolved across
    // the four inputs before opening large MC files
    const std::map<std::string, std::vector<std::string>> systDatasets =
        result.systematics.GetDataSets();
    const std::map<std::string, std::vector<std::string>> systParams =
        result.systematics.GetParamNames();
    for (std::map<std::string, std::vector<std::string>>::const_iterator entry =
             systDatasets.begin(); entry != systDatasets.end(); ++entry)
      for (const std::string &dataset : entry->second)
        if (!datasets.count(dataset))
          throw std::runtime_error("Systematic " + entry->first +
                                   " refers to unknown dataset " + dataset);
    for (std::map<std::string, std::vector<std::string>>::const_iterator entry =
             systParams.begin(); entry != systParams.end(); ++entry)
      for (const std::string &parameter : entry->second)
        if (!parameters.count(parameter))
          throw std::runtime_error("Systematic " + entry->first +
                                   " refers to unknown parameter " + parameter);

    return result;
  }
}
