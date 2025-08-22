#ifndef __ANTINUFIT__FitConfigLoader__
#define __ANTINUFIT__FitConfigLoader__

// Antinu headers
#include <FitConfig.hh>
#include <Utilities.hh>

// OXO headers
#include <ConfigLoader.hh>

// c++ headers
#include <algorithm>

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
