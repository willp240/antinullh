#pragma once

#include <string>
#include <stdexcept>

#include <StatisticSum.h>

namespace antinufit
{

  using ShapeFunc = std::function<double(const ParameterDict &)>;

  inline ShapeFunc GetShapeConstrFunc(const std::string &funcName)
  {
    if (funcName == "scaleFrac")
    {
      return [](const ParameterDict &params)
      {
        double total_rate_U =
            params.at("geonu_U_norm") *
            (params.at("geonu_U") + params.at("geonu_U2"));

        double total_rate_Th =
            params.at("geonu_Th_norm") *
            (params.at("geonu_Th") + params.at("geonu_Th2"));

        return (total_rate_U - total_rate_Th) /
               (total_rate_U + total_rate_Th);
      };
    }

    throw std::runtime_error("Unknown shape constraint function: " + funcName);
  }

} // namespace antinufit
