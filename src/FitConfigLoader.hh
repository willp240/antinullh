#ifndef __ANTINUFIT__FitConfigLoader__
#define __ANTINUFIT__FitConfigLoader__
#include <string>
#include <FitConfig.hh>
#include <map>

namespace antinufit
{
  class FitConfig;
  class FitConfigLoader
  {
  public:
    FitConfigLoader(const std::string &filePath_);
    ~FitConfigLoader();
    FitConfig LoadActive() const;

  private:
    std::string fPath;
  };
}
#endif
