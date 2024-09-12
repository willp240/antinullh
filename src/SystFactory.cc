#include <SystFactory.hh>

namespace antinufit
{

  Systematic *
  SystFactory::New(const std::string &name,
                   const std::string &type_,
                   const std::vector<std::string> &paramnamevec_,
                   ParameterDict &paramvals_,
                   std::string function_)
  {
    Systematic *syst;

    if (type_ == "scale")
    {
      Scale *scale = new Scale("scale");
      scale->RenameParameter("scaleFactor", paramnamevec_.at(0));
      scale->SetScaleFactor(paramvals_[paramnamevec_.at(0)]);
      syst = scale;
    }

    else if (type_ == "shift")
    {
      Shift *shift = new Shift("shift");
      shift->RenameParameter("shift", paramnamevec_.at(0));
      shift->SetShift(paramvals_[paramnamevec_.at(0)]);
      syst = shift;
    }

    else if (type_ == "sqroot_scale_conv")
    {
      VaryingCDF *smearer = new VaryingCDF("smear");
      // The 0 and 1.0 are arbitrary here. The parameter of the SquareRootScale is what will be a fit component
      Gaussian *gaus = new Gaussian(0, 1.0, paramnamevec_.at(0));
      smearer->SetKernel(gaus);

      SquareRootScale *sqrtscale = new SquareRootScale("e_smear_sigma_func");
      sqrtscale->RenameParameter("grad", paramnamevec_.at(0) );
      sqrtscale->SetGradient(paramvals_[paramnamevec_.at(0)]);
      smearer->SetDependance("stddevs_0", sqrtscale);

      Convolution *conv = new Convolution("conv");
      conv->SetConditionalPDF(smearer);
      syst = conv;
    }

    else if (type_ == "scale_function")
    {
      ScaleFunction *scale_func = new ScaleFunction("scale_function");
      auto func = std::get_if<std::function<double(const ParameterDict&, const double&)>>(&(functionMap[function_]));
      scale_func->SetScaleFunction(*func, paramnamevec_);
      scale_func->RenameParameter(paramnamevec_.at(0), "birks_constant");
      ParameterDict params({{"birks_constant", paramvals_[paramnamevec_.at(0)]}});
      scale_func->SetParameters(params);
      syst = scale_func;
    }

    else if (type_ == "shape")
    {
      Shape *shape = new Shape("shape");
      auto func = std::get_if<std::function<double(const ParameterDict&, const std::vector<double>&)>>(&(functionMap[function_]));
      shape->SetShapeFunction(*func, paramnamevec_);
      shape->RenameParameter(paramnamevec_.at(0), "deltam21");
      shape->RenameParameter(paramnamevec_.at(1), "theta12");
      ParameterDict params({{"deltam21", paramvals_[paramnamevec_.at(0)]}, {"theta12", paramvals_[paramnamevec_.at(1)]}});
      shape->SetParameters(params);
      syst = shape;
    }

    else
    {
      throw ValueError("Unknown systematic type: " + type_);
    }

    return syst;
  }
}
