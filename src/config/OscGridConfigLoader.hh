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
  class OscGridConfig;
  class OscGridConfigLoader
  {
  public:
    OscGridConfigLoader(const std::string &filePath_);
    ~OscGridConfigLoader();
    OscGridConfig Load() const;

  private:
    std::string fPath;
  };
}
#endif
