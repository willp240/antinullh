#pragma once

#include <cmath>
#include <stdexcept>
#include <string>

#include <ParameterDict.h>

namespace antinufit
{

  // Return which theta12-like parameter is present ("theta12", "sintheta12", or "sinsqtheta12").
  // Throws if none / more than one is present, or if deltam21 is missing.
  inline std::string GetTheta12ParamName(const ParameterDict &noms)
  {
    const bool hasTheta12 = noms.find("theta12") != noms.end();
    const bool hasSinTheta12 = noms.find("sintheta12") != noms.end();
    const bool hasSinSqTheta12 = noms.find("sinsqtheta12") != noms.end();
    const bool hasDeltam21 = noms.find("deltam21") != noms.end();

    const int thetacount = int(hasTheta12) + int(hasSinTheta12) + int(hasSinSqTheta12);

    if (!hasDeltam21 || thetacount == 0)
      throw std::runtime_error(
          "ERROR: A theta12, sintheta12, or sinsqtheta12 parameter, "
          "along with a deltam21 parameter, must be provided in the fitconfig.");

    if (thetacount > 1)
      throw std::runtime_error(
          "ERROR: More than one of theta12, sintheta12, sinsqtheta12 parameters "
          "were set in the fitconfig. There must be exactly one.");

    if (hasTheta12)
      return "theta12";
    if (hasSinTheta12)
      return "sintheta12";
    return "sinsqtheta12";
  }

  // Convert the stored theta parameter value to sin^2(theta12).
  // Assumes:
  //   - "theta12"      is in degrees
  //   - "sintheta12"   is sin(theta12)
  //   - "sinsqtheta12" is sin^2(theta12)
  inline double Theta12ToSinSq(const std::string &theta12name, double thetaValue)
  {
    if (theta12name == "sinsqtheta12")
      return thetaValue;

    if (theta12name == "sintheta12")
      return thetaValue * thetaValue;

    if (theta12name == "theta12")
    {
      const double thetaRad = M_PI * thetaValue / 180.0;
      const double s = std::sin(thetaRad);
      return s * s;
    }

    throw std::runtime_error("Unknown theta12 parameter name: " + theta12name);
  }

  // Bundle of oscillation-parameter conventions + nominal values.
  struct OscParams
  {
    std::string theta12name; // "theta12" | "sintheta12" | "sinsqtheta12"
    double deltam21 = 0.0;   // nominal (from noms)
    double theta12 = 0.0;    // nominal raw value (deg or sin or sinsq)

    double theta12_sinsq() const { return Theta12ToSinSq(theta12name, theta12); }
  };

  // Build OscParams from nominal values. Validates required params exist via GetTheta12ParamName.
  inline OscParams BuildOscParams(const ParameterDict &noms)
  {
    OscParams o;
    o.theta12name = GetTheta12ParamName(noms);
    o.deltam21 = noms.at("deltam21");
    o.theta12 = noms.at(o.theta12name);
    return o;
  }

} // namespace antinufit
