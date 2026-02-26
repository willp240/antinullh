#pragma once

// Antinu headers
#include <SystConfig.hh>

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
