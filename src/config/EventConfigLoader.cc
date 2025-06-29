#include <EventConfigLoader.hh>

namespace antinufit
{

  EventConfigLoader::EventConfigLoader(const std::string &filePath_)
  {
    fPath = filePath_;
  }

  EventConfig
  EventConfigLoader::LoadOne(const std::string &name_) const
  {

    ConfigLoader::Open(fPath);
    std::string baseDir;
    std::string prunedDir;

    ConfigLoader::Load("summary", "orig_base_dir", baseDir);
    ConfigLoader::Load("summary", "pruned_ntup_dir", prunedDir);

    std::vector<std::string> ntupFiles;
    int numDimensions;
    std::vector<std::string> groups;

    ConfigLoader::Load(name_, "ntup_files", ntupFiles);
    ConfigLoader::Load(name_, "dimensions", numDimensions);

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

  std::map<std::string, std::map<std::string, EventConfig> >
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
      if (std::find(toLoad.begin(), toLoad.end(), "all") != toLoad.end()){
        dsMap[*itDS] = LoadAll(dontLoad);
        continue;
      }

      EventConfigMap evMap;
      for (StringSet::iterator itEv = toLoad.begin(); itEv != toLoad.end(); ++itEv)
        evMap[*itEv] = LoadOne(*itEv);

      dsMap[*itDS] = evMap;
    }

    return dsMap;
  }

  EventConfigLoader::~EventConfigLoader()
  {
    ConfigLoader::Close();
  }

  std::map<std::string, EventConfig>
  EventConfigLoader::LoadAll(const std::set<std::string> &except_) const
  {
    typedef std::set<std::string> StringSet;
    StringSet toLoad = ConfigLoader::ListSections();
    toLoad.erase("summary");

    std::map<std::string, EventConfig> evMap;
    for (StringSet::iterator it = toLoad.begin(); it != toLoad.end();
         ++it)
    {
      if (!except_.count(*it))
        evMap[*it] = LoadOne(*it);
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
