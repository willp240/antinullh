#pragma once

#include <string>
#include <unordered_map>

#include <TVector3.h>

namespace antinufit
{
  TVector3 LLAtoECEF(double longitude, double latitude, double altitude);
  double GetReactorDistanceLLA(const double &longitude,
                               const double &latitude,
                               const double &altitude);
  std::unordered_map<int, double> LoadIndexDistanceMap(const std::string &filename);
  std::unordered_map<std::string, int> LoadNameIndexMap(const std::string &filename);
}
