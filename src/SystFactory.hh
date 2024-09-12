#ifndef __ANTINUFIT__SystFactory__
#define __ANTINUFIT__SystFactory__

#include <string>
#include <map>
#include <variant>

#include <Shift.h>
#include <Scale.h>
#include <Shape.h>
#include <SquareRootScale.h>
#include <Convolution.h>
#include <Gaussian.h>
#include <VaryingCDF.h>
#include <ScaleFunction.h>
#include <Exceptions.h>

#include <Functions.hh>

class Systematic;
namespace antinufit
{
  class SystFactory
  {
  public:
    static Systematic *New(const std::string &name,
                           const std::string &type_,
                           const std::vector<std::string> &paramnamevec_,
                           ParameterDict &paramvals_,
                           std::string function_ = "");
 };
}
#endif
