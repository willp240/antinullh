#ifndef __ANTINUFIT__Functions__
#define __ANTINUFIT__Functions__

// OXO headers
#include <ParameterDict.h>

// c++ headers
#include <vector>
#include <functional>
#include <variant>
#include <nlohmann/json.hpp>
#include <fstream>
#include <iostream>

// RAT headers
#include <RAT/DB.hh>

namespace antinufit
{
  double BirksLaw(const ParameterDict &, const double &);
  double OscProb(const ParameterDict &, const std::vector<double> &);
  // Disappearance probability 1 - P_ee(baseline [km], nuEnergy [MeV], deltam21, sin^2(theta12))
  double OscProb2(double, double, double, double );

  // OscProb2 split into its baseline-independent part (computed once per energy) and a cheap
  // per-baseline part: OscProb2EnergyTerms(E, dm21, s12sq).Evaluate(L) == OscProb2(L, E, dm21, s12sq)
  struct OscProb2Terms
  {
    double nuEnergy;
    double dEigen10, dEigen20, dEigen21;
    double X10, X20, X21;
    double Evaluate(double baseline) const;
  };
  OscProb2Terms OscProb2EnergyTerms(double nuEnergy, double fDmSqr21, double fSSqrTheta12);

  using FunctionVariant = std::variant<
      std::function<double(const ParameterDict &, const double &)>,
      std::function<double(const ParameterDict &, const std::vector<double> &)>>;

  extern std::map<std::string, FunctionVariant> functionMap;

}
#endif
