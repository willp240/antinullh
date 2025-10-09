#include <EventConfigLoader.hh>

namespace antinufit
{

  EventConfigLoader::EventConfigLoader(const std::string &filePath_)
  {
    fPath = filePath_;
  }

  EventConfig
  EventConfigLoader::LoadOne(const std::string &name_, const std::string dataset_) const
  {

    ConfigLoader::Open(fPath);

    std::vector<std::string> ntupFiles;
    int numDimensions;
    std::vector<std::string> groups;
    std::string baseDir;
    std::string prunedDir;

    ConfigLoader::Load(name_, "ntup_files", ntupFiles);
    ConfigLoader::Load(name_, "dimensions", numDimensions);

    ConfigLoader::Load(dataset_, "orig_base_dir", baseDir);
    ConfigLoader::Load(dataset_, "pruned_ntup_dir", prunedDir);

    try
    {
      ConfigLoader::Load(name_, "groups", groups);
    }
    catch (...)
    {
      groups.push_back("");
    }

    if (groups.size() == 0)
      groups.push_back("");

    EventConfig retVal;
    retVal.SetNtupFiles(ntupFiles);
    retVal.SetName(name_);
    retVal.SetNtupBaseDir(baseDir);
    retVal.SetPrunedPath(prunedDir + "/" + name_ + ".root");
    retVal.SetGroup(groups);
    retVal.SetNumDimensions(numDimensions);
    return retVal;
  }

  std::map<std::string, std::map<std::string, EventConfig>>
  EventConfigLoader::LoadActive() const
  {
    typedef std::set<std::string> StringSet;
    typedef std::map<std::string, std::map<std::string, EventConfig>> DataSetConfigMap;
    typedef std::map<std::string, EventConfig> EventConfigMap;
    ConfigLoader::Open(fPath);

    StringSet dataSets;
    DataSetConfigMap dsMap;
    ConfigLoader::Load("summary", "datasets", dataSets);

    for (StringSet::iterator itDS = dataSets.begin(); itDS != dataSets.end(); ++itDS)
    {

      StringSet toLoad;
      StringSet dontLoad;

      ConfigLoader::Load(*itDS, "active", toLoad);
      ConfigLoader::Load(*itDS, "inactive", dontLoad);

      //  if all is in the list, just do all of them
      if (std::find(toLoad.begin(), toLoad.end(), "all") != toLoad.end())
      {
        dsMap[*itDS] = LoadAll(dontLoad, *itDS);
        continue;
      }

      EventConfigMap evMap;
      for (StringSet::iterator itEv = toLoad.begin(); itEv != toLoad.end(); ++itEv)
        evMap[*itEv] = LoadOne(*itEv, *itDS);

      dsMap[*itDS] = evMap;
    }

    return dsMap;
  }

  std::map<std::string, std::map<std::string, EventConfig>>
  EventConfigLoader::LoadActiveAndData() const
  {
    typedef std::set<std::string> StringSet;
    std::map<std::string, std::map<std::string, EventConfig>> activeEvs = LoadActive();

    StringSet dataSets;
    std::string baseDir;
    std::string prunedDir;
    std::vector<std::string> dataFilename;

    ConfigLoader::Load("summary", "datasets", dataSets);

    for (StringSet::iterator itDS = dataSets.begin(); itDS != dataSets.end(); ++itDS)
    {

      std::string dataset = *itDS;

      ConfigLoader::Load(dataset, "data_file", dataFilename);
      ConfigLoader::Load(dataset, "orig_base_dir", baseDir);
      ConfigLoader::Load(dataset, "pruned_ntup_dir", prunedDir);

      EventConfig dataEveCfg;
      dataEveCfg.SetNtupFiles(dataFilename);
      dataEveCfg.SetName("data");
      dataEveCfg.SetNtupBaseDir(baseDir);
      dataEveCfg.SetPrunedPath(prunedDir + "/data.root");

      std::map<std::string, EventConfig> evMap;
      evMap["data"] = dataEveCfg;

      activeEvs[dataset] = evMap;
    }
  }

  EventConfigLoader::~EventConfigLoader()
  {
    ConfigLoader::Close();
  }

  std::map<std::string, EventConfig>
  EventConfigLoader::LoadAll(const std::set<std::string> &except_, const std::string dataset_) const
  {
    typedef std::set<std::string> StringSet;
    StringSet toLoad = ConfigLoader::ListSections();
    toLoad.erase("summary");

    std::map<std::string, EventConfig> evMap;
    for (StringSet::iterator evIt = toLoad.begin(); evIt != toLoad.end(); ++evIt)
    {
      if (!except_.count(*evIt))
      {
        try
        {
          std::string ntupFiles;
          ConfigLoader::Load(*evIt, "ntup_files", ntupFiles);
          evMap[*evIt] = LoadOne(*evIt, dataset_);
        }
        catch (const std::exception &e)
        {
          continue;
        }
      }
    }
    return evMap;
  }

  std::map<std::string, std::string>
  EventConfigLoader::GetDataPaths() const
  {
    std::map<std::string, std::string> dataPathMap;
    typedef std::set<std::string> StringSet;
    ConfigLoader::Open(fPath);

    StringSet dataSets;
    ConfigLoader::Load("summary", "datasets", dataSets);

    for (StringSet::iterator itDS = dataSets.begin(); itDS != dataSets.end(); ++itDS)
    {

      std::string dataPath;
      ConfigLoader::Load(*itDS, "datafile", dataPath);
      dataPathMap[*itDS] = dataPath;
    }

    return dataPathMap;
  }
}
