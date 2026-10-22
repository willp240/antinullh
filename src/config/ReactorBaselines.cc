#include <ReactorBaselines.hh>

#include <TMath.h>
#include <nlohmann/json.hpp>

#include <cmath>
#include <fstream>
#include <stdexcept>

namespace antinufit
{
  TVector3 LLAtoECEF(double longitude, double latitude, double altitude)
  {
    // WGS84 Earth model; coordinates are returned in kilometres
    const double toRad = TMath::Pi() / 180.0;
    const double earthRadius = 6378137.0;
    const double flattening = 1.0 / 298.257223563;
    const double reducedLatitude =
        std::atan(std::pow(1.0 - flattening, 2) * std::tan(latitude * toRad));
    const double radius = std::sqrt(
        earthRadius * earthRadius /
        (1.0 + (1.0 / std::pow(1.0 - flattening, 2) - 1.0) *
                   std::pow(std::sin(reducedLatitude), 2)));
    const double lon = longitude * toRad;
    const double lat = latitude * toRad;
    const double x = (radius * std::cos(reducedLatitude) * std::cos(lon) +
                      altitude * std::cos(lat) * std::cos(lon)) / 1000.0;
    const double y = (radius * std::cos(reducedLatitude) * std::sin(lon) +
                      altitude * std::cos(lat) * std::sin(lon)) / 1000.0;
    const double z = (radius * std::sin(reducedLatitude) +
                      altitude * std::sin(lat)) / 1000.0;
    return TVector3(x, y, z);
  }

  double GetReactorDistanceLLA(const double &longitude,
                               const double &latitude,
                               const double &altitude)
  {
    const TVector3 sno(672.87, -4347.18, 4600.51);
    return (LLAtoECEF(longitude, latitude, altitude) - sno).Mag();
  }

  namespace
  {
    nlohmann::json ReadReactorJson(const std::string &filename)
    {
      std::ifstream file(filename);
      if (!file)
        throw std::runtime_error("Cannot open reactor JSON file: " + filename);
      nlohmann::json data;
      file >> data;
      return data;
    }
  }

  std::unordered_map<int, double> LoadIndexDistanceMap(const std::string &filename)
  {
    const nlohmann::json data = ReadReactorJson(filename);
    std::unordered_map<int, double> distances;
    for (nlohmann::json::const_iterator entry = data.cbegin(); entry != data.cend(); ++entry)
      distances.emplace(std::stoi(entry.key()), entry.value().at(1).get<double>());
    return distances;
  }

  std::unordered_map<std::string, int> LoadNameIndexMap(const std::string &filename)
  {
    const nlohmann::json data = ReadReactorJson(filename);
    std::unordered_map<std::string, int> indices;
    for (nlohmann::json::const_iterator entry = data.cbegin(); entry != data.cend(); ++entry)
      indices.emplace(entry.value().at(0).get<std::string>(), std::stoi(entry.key()));
    return indices;
  }
}
