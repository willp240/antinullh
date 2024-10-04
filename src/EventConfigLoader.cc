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
    std::string pdfDir;

    ConfigLoader::Load("summary", "orig_base_dir", baseDir);
    ConfigLoader::Load("summary", "pruned_ntup_dir", prunedDir);
    ConfigLoader::Load("summary", "pdf_dir", pdfDir);

    double rate;
    std::string texLabel;
    std::vector<std::string> ntupFiles;
    int numDimensions;
    std::string group;

    ConfigLoader::Load(name_, "rate", rate);
    ConfigLoader::Load(name_, "tex_label", texLabel);
    ConfigLoader::Load(name_, "ntup_files", ntupFiles);
    ConfigLoader::Load(name_, "dimensions", numDimensions);

    try
    {
      ConfigLoader::Load(name_, "groups", group);
    }
    catch (...)
    {
      group = "";
    }

    std::vector<std::string> groups;
    std::stringstream ss(group);
    std::string indiv_group;
    while (std::getline(ss, indiv_group, ','))
      groups.push_back(indiv_group);

    EventConfig retVal;
    retVal.SetRate(rate);
    retVal.SetNtupFiles(ntupFiles);
    retVal.SetTexLabel(texLabel);
    retVal.SetName(name_);
    retVal.SetNtupBaseDir(baseDir);
    retVal.SetPrunedPath(prunedDir + "/" + name_ + ".root");
    retVal.SetPdfPath(baseDir + "/" + name_ + ".root");
    retVal.SetGroup(groups);
    retVal.SetNumDimensions(numDimensions);
    return retVal;
  }

  std::map<std::string, EventConfig>
  EventConfigLoader::LoadActive() const
  {
    typedef std::set<std::string> StringSet;
    typedef std::map<std::string, EventConfig> EventConfigMap;
    ConfigLoader::Open(fPath);

    StringSet toLoad;
    StringSet dontLoad;

    ConfigLoader::Load("summary", "active", toLoad);
    ConfigLoader::Load("summary", "inactive", dontLoad);
    //  if all is in the list, just do all of them
    if (std::find(toLoad.begin(), toLoad.end(), "all") != toLoad.end())
      return LoadAll(dontLoad);

    EventConfigMap evMap;
    for (StringSet::iterator it = toLoad.begin(); it != toLoad.end(); ++it)
      evMap[*it] = LoadOne(*it);

    return evMap;
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

}
