#pragma once

// Antinu headers
#include <Functions.hh>
#include <Utilities.hh>

// OXO headers
#include <Shift.h>
#include <Scale.h>
#include <Shape.h>
#include <SquareRootScale.h>
#include <Convolution.h>
#include <Gaussian.h>
#include <VaryingCDF.h>
#include <ScaleFunction.h>
#include <Exceptions.h>

// c++ headers
#include <string>
#include <map>
#include <variant>

class Systematic;
namespace antinufit
{
  class SystFactory
  {
  public:
    static Systematic *New(const std::string &name,
                           const std::string &type_,
                           const std::vector<std::string> &paramnamevec_,
                           const ParameterDict &paramvals_,
                           const ShapeFunction &shapeFunction_ = ShapeFunction());
  };

}
