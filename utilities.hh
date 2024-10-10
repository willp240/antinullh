#ifndef __ANTINUFIT__Utilities__
#define __ANTINUFIT__Utilities__

// c++ headers
#include <vector>
#include <fstream>
#include <nlohmann/json.hpp>

namespace antinufit
{

  std::vector<double> linspace(double, double, size_t);
  std::unordered_map<int, double> LoadIndexDistanceMap(std::string );
}
#endif
