#pragma once

// Antinu headers
#include <OscGridConfig.hh>

// OXO headers
#include <ConfigLoader.hh>

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
