#ifndef __ANTINUFIT__DistConfigLoader__
#define __ANTINUFIT__DistConfigLoader__
#include <string>
#include <DistConfig.hh>
#include <map>

namespace antinufit{
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
