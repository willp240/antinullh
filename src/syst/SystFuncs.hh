#pragma once

// OXO headers
#include <ParameterDict.h>
#include <BinnedED.h>
#include <Shape.h>

// c++ headers
#include <vector>
#include <unordered_map>
#include <nlohmann/json.hpp>
#include <fstream>
#include <iostream>

// RAT headers
#include <RAT/DB.hh>

namespace antinufit
{
  double BirksLaw(double energy, double birksConstant);
  ShapeFunction MakeAlphaNClassifierShape(const std::string &functionName,
                                          const std::string &parameterName);
  // Disappearance probability 1 - P_ee(baseline [km], nuEnergy [MeV], deltam21, sin^2(theta12))
  double OscProb(double, double, double, double);

  // Prepare a Shape function from MC binned in (energy, reactorIndex, nu_energy)
  // It returns oscillated/original counts for the requested energy bin; empty bins have weight 1
  ShapeFunction MakeReactorOscillationShape(const BinnedED &reactorMC,
                                           const std::unordered_map<int, double> &baselines,
                                           const std::string &theta12Name);
}
