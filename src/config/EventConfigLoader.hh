#ifndef __ANTINUFIT__EventConfigLoader__
#define __ANTINUFIT__EventConfigLoader__

// Antinu headers
#include <EventConfig.hh>

// OXO headers
#include <ConfigLoader.hh>

// c++ headers
#include <algorithm>

namespace antinufit
{
  class EventConfig;
  class EventConfigLoader
  {
  public:
    EventConfigLoader(const std::string &filePath_);
    ~EventConfigLoader();
    EventConfig LoadOne(const std::string &name_) const;
    std::map<std::string, std::map<std::string, EventConfig>> LoadActive() const;
    std::map<std::string, EventConfig> LoadAll(const std::set<std::string> &except_) const;
    std::map<std::string, std::string> GetDataPaths() const;

  private:
    std::string fPath;
  };
}
#endif
