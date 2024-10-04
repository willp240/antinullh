#include <SystConfigLoader.hh>

namespace antinufit
{

  SystConfigLoader::SystConfigLoader(const std::string &filePath_)
  {
    fPath = filePath_;
  }

  SystConfig
  SystConfigLoader::LoadActive() const
  {
    SystConfig ret;
    ConfigLoader::Open(fPath);

    typedef std::set<std::string> StringSet;
    StringSet toLoad;
    ConfigLoader::Load("summary", "active", toLoad);

    std::string name;
    std::string paramNames;
    std::vector<std::string> distObs;
    std::vector<std::string> transObs;
    std::string type;
    std::string group;
    std::string function;

    if (std::find(toLoad.begin(), toLoad.end(), "all") != toLoad.end())
    {
      toLoad = ConfigLoader::ListSections();
      toLoad.erase("summary");
    }

    for (StringSet::iterator it = toLoad.begin(); it != toLoad.end(); ++it)
    {
      name = *it;
      ConfigLoader::Load(name, "dist_obs", distObs);
      ConfigLoader::Load(name, "trans_obs", transObs);
      ConfigLoader::Load(name, "type", type);
      try
      {
        ConfigLoader::Load(name, "group", group);
      }
      catch (ConfigFieldMissing)
      {
        group = "";
      }
      try
      {
        ConfigLoader::Load(name, "param_names", paramNames);
      }
      catch (ConfigFieldMissing)
      {
        paramNames = name;
      }
      if(type == "scale_function" || type == "shape"){
          try
          {
            ConfigLoader::Load(name, "function", function);
          }
          catch (ConfigFieldMissing) //FIXME change this oor
          {
            throw NotFoundError(Formatter() << "SystConfigLoader::No function set for " << type << " systematic: " << name);
          }
          if(functionMap.find(function) == functionMap.end())
            throw NotFoundError(Formatter() << "SystConfigLoader::Invalid function ( " << function << ") set for systematic: " << name);
      }
      else
        function == "";

      ret.AddParameter(name, paramNames, distObs, transObs, type, group, function);
      distObs.clear();
      transObs.clear();
    }
    return ret;
  }

  SystConfigLoader::~SystConfigLoader()
  {
    ConfigLoader::Close();
  }
}
