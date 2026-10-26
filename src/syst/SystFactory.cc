#include <SystFactory.hh>
#include <OscParams.hh>
#include <set>

namespace antinufit
{

  Systematic *
  SystFactory::New(const std::string &name,
                   const std::string &type_,
                   const std::vector<std::string> &paramnamevec_,
                   const ParameterDict &paramvals_,
                   const ShapeFunction &shapeFunction_)
  {
    Systematic *syst;

    std::vector<std::string> type_vec = SplitString(type_, ':');
    std::string type = type_vec.size() > 0 ? type_vec[0] : "";
    std::string function = type_vec.size() > 1 ? type_vec[1] : "";
    std::string ploy = type_vec.size() > 2 ? type_vec[2] : "";

    if (type == "Scale")
    {
      Scale *scale = new Scale(name);
      scale->RenameParameter("scaleFactor", paramnamevec_.at(0));
      scale->SetScaleFactor(paramvals_.at(paramnamevec_.at(0)));
      syst = scale;
    }

    else if (type == "Shift")
    {
      Shift *shift = new Shift(name);
      shift->RenameParameter("shift", paramnamevec_.at(0));
      shift->SetShift(paramvals_.at(paramnamevec_.at(0)));
      syst = shift;
    }

    else if (type == "Conv")
    {
      // First declare possible functions and ploys
      // The 0 and 1.0 are arbitrary here. The parameter of the SquareRootScale is what will be a fit component
      Gaussian *gaus = new Gaussian(0, 1.0, "gaus");
      SquareRootScale *sqrtscale = new SquareRootScale("e_smear_sigma_func");

      VaryingCDF *smearer = new VaryingCDF("smear");

      if (function == "Gaussian")
      {
        smearer->SetKernel(gaus);
      }
      else
      {
        throw ValueError("Unknown function, " + function + ", for systematic type: " + type);
      }

      if (ploy == "SquareRootScale")
      {
        sqrtscale->RenameParameter("grad", paramnamevec_.at(0));
        sqrtscale->SetGradient(paramvals_.at(paramnamevec_.at(0)));
        smearer->SetDependance("stddevs_0", sqrtscale);
      }
      else
      {
        throw ValueError("Unknown ploy, " + ploy + ", for systematic type: " + type);
      }

      Convolution *conv = new Convolution(name);
      conv->SetConditionalPDF(smearer);
      syst = conv;
    }

    else if (type == "ScaleFunction")
    {
      ScaleFunction *scale_func = new ScaleFunction(name);
      if (function == "BirksLaw" || function == "BirksLaw2")
      {
        const std::string parameter = paramnamevec_.at(0);
        const ScaleFunc birks = [parameter](const ParameterDict &params, const double &energy)
        {
          return BirksLaw(energy, params.at(parameter));
        };
        scale_func->SetScaleFunction(birks, paramnamevec_);
        ParameterDict params({{parameter, paramvals_.at(paramnamevec_.at(0))}});
        scale_func->SetParameters(params);
      }
      else
      {
        throw ValueError("Unknown function, " + function + ", for systematic type: " + type);
      }
      syst = scale_func;
    }

    else if (type == "Shape")
    {

      static const std::set<std::string> normalisationNames = {
          "reactor_nubar_norm", "geonu_U_norm", "geonu_Th_norm",
          "alphan_PRecoil_norm", "alphan_CScatter_norm", "alphan_OExcited_norm",
          "bipolike_norm", "atmospheric_norm"};
      Shape *shape = new Shape(name);

      if (shapeFunction_)
      {
        shape->SetShapeFunction(shapeFunction_, paramnamevec_);
        for (const std::string &param : paramnamevec_)
          shape->SetParameter(param, paramvals_.at(param));
      }
      else if (function == "GeoOscProb")
      {
        shape->SetShapeFunction(MakeGeoOscProbShape(paramnamevec_), paramnamevec_);
        shape->RenameParameter(paramnamevec_.at(0), "deltam21");
        shape->RenameParameter(paramnamevec_.at(1), "sinsqtheta12");
        ParameterDict params({{"deltam21", paramvals_.at(paramnamevec_.at(0))}, {"sinsqtheta12", paramvals_.at(paramnamevec_.at(1))}});
        shape->SetParameters(params);
      }

      else if (function.rfind("AlphaNClass", 0) == 0)
      {
        const bool alphaN = function == "AlphaNClassAlphaNPPO" ||
                            function == "AlphaNClassAlphaNBisMSB";
        const bool bisMSB = function.find("BisMSB") != std::string::npos;
        const std::string parameter = (alphaN ? "class_a_" : "class_s_") +
                                      std::string(bisMSB ? "bismsb" : "ppo");
        shape->SetShapeFunction(MakeAlphaNClassifierShape(function, parameter), paramnamevec_);
        shape->RenameParameter(paramnamevec_.at(0), parameter);
        ParameterDict params({{parameter, paramvals_.at(paramnamevec_.at(0))}});
        shape->SetParameters(params);
      }

      else if (normalisationNames.count(function))
      {
        const ShapeFunction normalisation = [function](const ParameterDict &params,
                                                       const std::vector<double> &)
        {
          return params.at(function);
        };
        shape->SetShapeFunction(normalisation, paramnamevec_);
        shape->RenameParameter(paramnamevec_.at(0), function);
        ParameterDict params({{function, paramvals_.at(paramnamevec_.at(0))}});
        shape->SetParameters(params);
      }
      else
      {
        throw ValueError("Unknown function, " + function + ", for systematic type: " + type);
      }

      syst = shape;
    }

    else
    {
      throw ValueError("Unknown systematic type: " + type_);
    }

    return syst;
  }
}
