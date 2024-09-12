#ifndef __ANTINUFIT__Functions__
#define __ANTINUFIT__Functions__

#include <variant>
#include <ParameterDict.h>
#include <vector>
#include <functional>

namespace antinufit
{
  double BirksLaw(const ParameterDict&, const double&);
  double OscProb(const ParameterDict&, const std::vector<double>&);

  using FunctionVariant = std::variant<
      std::function<double(const ParameterDict&, const double&)>,
      std::function<double(const ParameterDict&, const std::vector<double>&)>>;

  extern std::map<std::string, FunctionVariant> functionMap;

}
#endif
