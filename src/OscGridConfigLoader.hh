#ifndef __ANTINUFIT__OscGridConfigLoader__
#define __ANTINUFIT__OscGridConfigLoader__

// Antinu headers
#include <OscGridConfig.hh>

// OXO headers
#include <ConfigLoader.hh>

// c++ headers
#include <algorithm>

namespace antinufit
{
  class OscConfig;
  class OscConfigLoader
  {
  public:
    OscConfigLoader(const std::string &filePath_);
    ~FitConfigLoader();

  private:
    std::string fPath;
  };
}
#endif
