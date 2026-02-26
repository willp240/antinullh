#pragma once

#include <map>
#include <string>
#include <vector>

#include <EventConfig.hh>
#include <PDFConfig.hh>
#include <EventConfigLoader.hh>
#include <PDFConfigLoader.hh>
#include <DistBuilder.hh>

namespace antinufit
{

  using EvMap = std::map<std::string, EventConfig>;
  using DSMap = std::map<std::string, EvMap>;

  // Info for building the event dists for each dataset
  struct DatasetSetup
  {
    DSMap dsPDFMap;
    std::map<std::string, std::string> dataPath;

    PDFConfig pdfConfig;
    std::vector<std::string> dataObs;
    ObsSet dataObsSet;
    AxisCollection systAxes;
  };

  inline DatasetSetup LoadDatasetSetup(const std::string &evConfigFile,
                                       const std::string &pdfConfigFile)
  {
    DatasetSetup setup;

    // Event config
    EventConfigLoader evLoader(evConfigFile);
    setup.dsPDFMap = evLoader.LoadActive();
    setup.dataPath = evLoader.GetPrunedDataPaths();

    // PDF config
    PDFConfigLoader pdfLoader(pdfConfigFile);
    setup.pdfConfig = pdfLoader.Load();
    setup.dataObs = setup.pdfConfig.GetDataBranchNames();
    setup.dataObsSet = ObsSet(setup.dataObs);
    setup.systAxes = DistBuilder::BuildAxes(setup.pdfConfig,
                                            setup.dataObs.size());

    return setup;
  }

} // namespace antinufit
