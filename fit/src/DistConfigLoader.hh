#ifndef __BBFIT__DistConfigLoader__
#define __BBFIT__DistConfigLoader__
#include <string>
#include <DistConfig.hh>
#include <map>

namespace bbfit{
class DistConfig;
class DistConfigLoader{
public:
  DistConfigLoader(const std::string& filePath_);
  ~DistConfigLoader();
  DistConfig Load() const;

private:
  std::string fPath;
};
}
#endif
