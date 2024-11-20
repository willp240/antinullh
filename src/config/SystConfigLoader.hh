#ifndef __ANTINUFIT__SystConfigLoader__
#define __ANTINUFIT__SystConfigLoader__

// Antinu headers
#include <SystConfig.hh>
#include <Functions.hh>

// OXO headers
#include <ConfigLoader.hh>

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
