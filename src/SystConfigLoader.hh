#ifndef __ANTINUFIT__SystConfigLoader__
#define __ANTINUFIT__SystConfigLoader__
#include <string>
#include <SystConfig.hh>
#include <vector>
#include <Functions.hh>

namespace antinufit
{
  class SystConfig;
  class SystConfigLoader
  {
  public:
    SystConfigLoader(const std::string &filePath_);
    ~SystConfigLoader();
    SystConfig LoadActive() const;

  private:
    std::string fPath;
  };
}
#endif
