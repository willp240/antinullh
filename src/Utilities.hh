#ifndef __ANTINUFIT__Utilities__
#define __ANTINUFIT__Utilities__

// C++ headers
#include <vector>
#include <fstream>
#include <nlohmann/json.hpp>
#include <iostream>

// ROOT headers
#include <TVector3.h>

namespace antinufit
{

  TVector3 LLAtoECEF(double, double, double);
  double GetReactorDistanceLLA(const double &, const double &, const double &);
  std::unordered_map<int, double> LoadIndexDistanceMap(std::string);
  std::unordered_map<std::string, int> LoadNameIndexMap(std::string);
  std::vector<double> linspace(double, double, size_t);
  std::pair<size_t, size_t> GetLowerUpperIndices(const std::vector<double>, double);
}
#endif
