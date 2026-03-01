#pragma once

// Antinu headers
#include <FitConfig.hh>

// OXO headers
#include <ConfigLoader.hh>

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
