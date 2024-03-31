#ifndef __BBFIT__CutConfigLoader__
#define __BBFIT__CutConfigLoader__
#include <string>
#include <CutConfig.hh>
#include <vector>

namespace bbfit{
class CutConfig;
class CutConfigLoader{
public:
  CutConfigLoader(const std::string& filePath_);
  ~CutConfigLoader();
  CutConfig LoadOne(const std::string& name_) const;
  std::vector<CutConfig> LoadActive() const;

private:
  std::string fPath;
};
}
#endif
