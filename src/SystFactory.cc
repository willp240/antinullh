#include <SystFactory.hh>
#include <Shift.h>
#include <Scale.h>
#include <SquareRootScale.h>
#include <Convolution.h>
#include <Gaussian.h>
#include <VaryingCDF.h>
#include <ScaleFunction.h>
#include <Exceptions.h>

double birks_law(const ParameterDict &params, const double &obs_val)
{
  const double birks_const = 0.074;
  return ((1 + (birks_const * obs_val)) / (1 + (params.at("birks_constant") * obs_val))) * obs_val;
}

namespace antinufit
{

  Systematic *
  SystFactory::New(const std::string &name_,
                   const std::string &type_,
                   double value)
  {
    Systematic *syst;

    if (type_ == "scale")
    {
      Scale *scale = new Scale("scale");
      scale->RenameParameter("scaleFactor", name_);
      scale->SetScaleFactor(value);
      syst = scale;
    }

    else if (type_ == "shift")
    {
      Shift *shift = new Shift("shift");
      shift->RenameParameter("shift", name_);
      shift->SetShift(value);
      syst = shift;
    }

    else if (type_ == "sqroot_scale_conv")
    {
      VaryingCDF *smearer = new VaryingCDF("smear");
      // The 0 and 1.0 are arbitrary here. The parameter of the SquareRootScale is what will be a fit component
      Gaussian *gaus = new Gaussian(0, 1.0, name_);
      smearer->SetKernel(gaus);

      SquareRootScale *sqrtscale = new SquareRootScale("e_smear_sigma_func");
      sqrtscale->RenameParameter("grad", name_ + "_mean");
      sqrtscale->SetGradient(value);
      smearer->SetDependance("stddevs_0", sqrtscale);

      Convolution *conv = new Convolution("conv");
      conv->SetConditionalPDF(smearer);
      syst = conv;
    }

    else if (type_ == "scale_function")
    {
      ScaleFunction *scale_func = new ScaleFunction("scale_function");
      const std::vector<std::string> scale_param_names{name_};
      scale_func->SetScaleFunction(birks_law, scale_param_names);
      scale_func->RenameParameter(name_, "birks_constant");
      ParameterDict params({{"birks_constant", value}});
      scale_func->SetParameters(params);
      syst = scale_func;
    }

    else
    {
      throw ValueError("Unknown systematic type: " + type_);
    }

    return syst;
  }

}
