#ifndef __BBFIT__EventConfigLoader__
#define __BBFIT__EventConfigLoader__
#include <EventConfig.hh>
#include <map>
#include <string>
#include <set>

namespace bbfit{
class EventConfig;
class EventConfigLoader{
public:
  EventConfigLoader(const std::string& filePath_);
  ~EventConfigLoader();
  EventConfig LoadOne(const std::string& name_) const;
  std::map<std::string, EventConfig> LoadActive() const;
  std::map<std::string, EventConfig> LoadAll(const std::set<std::string>& except_) const;

private:
  std::string fPath;
};
}
#endif
