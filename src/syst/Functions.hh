#pragma once

// OXO headers
#include <ParameterDict.h>

// c++ headers
#include <vector>
#include <nlohmann/json.hpp>
#include <fstream>
#include <iostream>

// RAT headers
#include <RAT/DB.hh>

namespace antinufit
{
  double BirksLaw(const ParameterDict &, const double &);
  // Disappearance probability 1 - P_ee(baseline [km], nuEnergy [MeV], deltam21, sin^2(theta12))
  double OscProb(double, double, double, double);

  // OscProb split into its baseline-independent part (computed once per energy) and a cheap
  // per-baseline part: OscProbEnergyTerms(E, dm21, s12sq).Evaluate(L) is numerically equivalent
  // to OscProb(L, E, dm21, s12sq), up to floating-point rounding.
  struct OscProbTerms
  {
    double nuEnergy;
    double dEigen10, dEigen20, dEigen21;
    double X10, X20, X21;
    double Evaluate(double baseline) const;
  };
  OscProbTerms OscProbEnergyTerms(double nuEnergy, double fDmSqr21, double fSSqrTheta12);
}
