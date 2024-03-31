#ifndef __BBFIT__FitConfigLoader__
#define __BBFIT__FitConfigLoader__
#include <string>
#include <FitConfig.hh>
#include <map>

namespace bbfit{
class FitConfig;
class FitConfigLoader{
public:
  FitConfigLoader(const std::string& filePath_);
  ~FitConfigLoader();
  FitConfig LoadActive() const;
  
private:
  std::string fPath;
};
}
#endif
