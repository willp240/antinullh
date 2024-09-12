#include "Functions.hh"

namespace antinufit
{

  double BirksLaw(const ParameterDict &params, const double &obs_val)
  {
    const double birks_const = 0.074;
    return ((1 + (birks_const * obs_val)) / (1 + (params.at("birks_constant") * obs_val))) * obs_val;
  }

  double OscProb(const ParameterDict &params, const std::vector<double> &obs_vals)
  {
    // Calc osc prob as a function of params["theta12"], params["delatm21"], obs_val (energy)
    return 0;
  }

  std::map<std::string, FunctionVariant> functionMap = {
      {"BirksLaw", BirksLaw},
      {"OscProb", OscProb}
  };

}
