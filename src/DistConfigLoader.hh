#ifndef __ANTINUFIT__DistConfigLoader__
#define __ANTINUFIT__DistConfigLoader__

// Antinu headers
#include <DistConfig.hh>

// OXO headers
#include <ConfigLoader.hh>

// c++ headers
#include <algorithm>

namespace antinufit
{
  class DistConfig;
  class DistConfigLoader
  {
  public:
    DistConfigLoader(const std::string &filePath_);
    ~DistConfigLoader();
    DistConfig Load() const;

  private:
    std::string fPath;
  };
}
#endif
