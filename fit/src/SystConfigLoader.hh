#ifndef __BBFIT__SystConfigLoader__
#define __BBFIT__SystConfigLoader__
#include <string>
#include <SystConfig.hh>
#include <vector>

namespace bbfit{
class SystConfig;
class SystConfigLoader{
public:
  SystConfigLoader(const std::string& filePath_);
  ~SystConfigLoader();
  SystConfig LoadActive() const;

private:
  std::string fPath;
};
}
#endif
