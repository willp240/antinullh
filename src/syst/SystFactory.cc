#include <SystFactory.hh>

namespace antinufit
{

  Systematic *
  SystFactory::New(const std::string &name,
                   const std::string &type_,
                   const std::vector<std::string> &paramnamevec_,
                   ParameterDict &paramvals_)
  {
    std::map<int, OscGrid *> defaultGridMap;
    std::unordered_map<int, double> defaultIndexMap;
    return New(name, type_, paramnamevec_, paramvals_, defaultGridMap, defaultIndexMap);
  }
  Systematic *
  SystFactory::New(const std::string &name,
                   const std::string &type_,
                   const std::vector<std::string> &paramnamevec_,
                   ParameterDict &paramvals_,
                   std::map<int, OscGrid *> &oscgridmap_,
                   std::unordered_map<int, double> &indexdistancemap_)
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
      scale->SetScaleFactor(paramvals_[paramnamevec_.at(0)]);
      syst = scale;
    }

    else if (type == "Shift")
    {
      Shift *shift = new Shift(name);
      shift->RenameParameter("shift", paramnamevec_.at(0));
      shift->SetShift(paramvals_[paramnamevec_.at(0)]);
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
        sqrtscale->SetGradient(paramvals_[paramnamevec_.at(0)]);
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

      // First declare possible functions
      ScaleFunc BirksLaw = [](const ParameterDict &params, const double &obs_val)
      {
        const double birks_const = 0.074;
        return ((1 + (birks_const * obs_val)) / (1 + (params.at("birks_constant") * obs_val))) * obs_val;
      };

      ScaleFunc BirksLaw2 = [](const ParameterDict &params, const double &obs_val)
      {
        const double birks_const = 0.074;
        return ((1 + (birks_const * obs_val)) / (1 + (params.at("birks_constant2") * obs_val))) * obs_val;
      };

      ScaleFunction *scale_func = new ScaleFunction(name);

      if (function == "BirksLaw")
      {

        scale_func->SetScaleFunction(BirksLaw, paramnamevec_);
        scale_func->RenameParameter(paramnamevec_.at(0), "birks_constant");
        ParameterDict params({{"birks_constant", paramvals_[paramnamevec_.at(0)]}});
        scale_func->SetParameters(params);
      }
      else if (function == "BirksLaw2")
      {

        scale_func->SetScaleFunction(BirksLaw2, paramnamevec_);
        scale_func->RenameParameter(paramnamevec_.at(0), "birks_constant2");
        ParameterDict params({{"birks_constant2", paramvals_[paramnamevec_.at(0)]}});
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

      // First declare possible functions
      ShapeFunction OscProbGrid = [&oscgridmap_, &indexdistancemap_](const ParameterDict &params, const std::vector<double> &obs_vals)
      {
        double distance = indexdistancemap_[obs_vals.at(1)];
        double nuEnergy = obs_vals.at(2);
        OscGrid *oscGrid = oscgridmap_[obs_vals.at(1)];
        double prob = oscGrid->Evaluate(nuEnergy, params.at("deltam21"), params.at("theta12"));
        return prob;
      };

      ShapeFunction OscProb = [&indexdistancemap_](const ParameterDict &params, const std::vector<double> &obs_vals)
      {
        Double_t nuE_parent = obs_vals.at(2);
        Double_t baseline = indexdistancemap_[obs_vals.at(1)];
        Double_t fDmSqr21 = params.at("deltam21");
        Double_t fDmSqr32 = 2.451e-3;
        Double_t fSSqrTheta12 = sin(params.at("theta12")) * sin(params.at("theta12"));
        Double_t fSSqrTheta13 = 0.0216;

        // Declare quantities to use in loops
        Double_t scale;
        Double_t fOscProb = 0.0;
        Double_t Ne = 8.13e23;

        Double_t fDmSqr31 = fDmSqr32 + fDmSqr21;

        Double_t H_ee_vac = fDmSqr21 * (fSSqrTheta12 * (1 - fSSqrTheta13) - (1.0 / 3.0)) + fDmSqr31 * (fSSqrTheta13 - (1.0 / 3.0));
        Double_t H_neq2 = (1 - fSSqrTheta13) * (fDmSqr21 * fDmSqr21 * fSSqrTheta12 * (1 + fSSqrTheta12 * (fSSqrTheta13 - 1)) + fDmSqr31 * fDmSqr31 * fSSqrTheta13 - 2.0 * fDmSqr21 * fDmSqr31 * fSSqrTheta12 * fSSqrTheta13);

        Double_t a0_vac = -(2.0 / 27.0) * (fDmSqr21 * fDmSqr21 * fDmSqr21 + fDmSqr31 * fDmSqr31 * fDmSqr31) + (1.0 / 9.0) * (fDmSqr21 * fDmSqr21 * fDmSqr31 + fDmSqr21 * fDmSqr31 * fDmSqr31);
        Double_t a1_vac = (1.0 / 3.0) * (fDmSqr21 * fDmSqr31 - fDmSqr21 * fDmSqr21 - fDmSqr31 * fDmSqr31);
        Double_t Y_ee_vac = (2.0 / 3.0) * a1_vac + H_ee_vac * H_ee_vac + H_neq2;

        // Declare quantities to use in loop
        Double_t alpha = -2.535e-31 * Ne; // conversion factor in eV2/MeV
        Double_t A_CC;
        Double_t alpha_1;
        Double_t a0;
        Double_t a1;
        Double_t Y_ee;
        Double_t H_ee;
        Double_t eigen[3];
        Double_t X[3];
        Double_t arcCos;
        Double_t preFact;
        Double_t s_10;
        Double_t s_20;
        Double_t s_21;

        // do calculation
        scale = 1.267e3 * baseline / nuE_parent; // for nuE in [MeV] and baseline in [km]
        A_CC = alpha * nuE_parent;               // for A_CC in [eV^2] and nuE in [MeV]

        // Compute new values for H_ee, Y, a0 and a1 (make sure and Y are updated after their use by others)
        alpha_1 = H_ee_vac * A_CC + (1.0 / 3.0) * A_CC * A_CC;

        a0 = a0_vac - Y_ee_vac * A_CC - (1.0 / 3.0) * H_ee_vac * A_CC * A_CC - (2.0 / 27.0) * A_CC * A_CC * A_CC;
        a1 = a1_vac - alpha_1;
        Y_ee = Y_ee_vac + (2.0 / 3.0) * alpha_1;
        H_ee = H_ee_vac + (2.0 / 3.0) * A_CC;

        // Get eigenvalues of H, and constants X and theta
        arcCos = (1.0 / 3.0) * acos(1.5 * (a0 / a1) * sqrt(-3.0 / a1));
        preFact = 2.0 * sqrt(-a1 / 3.0);

        for (int i = 0; i < 3; ++i)
        {
          eigen[i] = preFact * cos(arcCos - (2.0 / 3.0) * M_PI * i);
          X[i] = (1.0 / 3.0) + (eigen[i] * H_ee + Y_ee) / (3.0 * eigen[i] * eigen[i] + a1);
        }

        s_10 = sin(scale * (eigen[1] - eigen[0]));
        s_20 = sin(scale * (eigen[2] - eigen[0]));
        s_21 = sin(scale * (eigen[2] - eigen[1]));

        // Compute probability
        fOscProb = 4.0 * (X[1] * X[0] * s_10 * s_10 + X[2] * X[0] * s_20 * s_20 + X[2] * X[1] * s_21 * s_21);

        return fOscProb;
      };

      // The bin-by-bin efficiencies scaling for the alpha,n classifier on the reactor IBD PDF
      ShapeFunction AlphaNClassIBD = [](const ParameterDict &params, const std::vector<double> &obs_vals)
      {
        // Nominal efficiency scaling for each bin
        std::vector<double> efficiencies = {0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8,
                                            0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8,
                                            0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8,
                                            0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8,
                                            0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8};

        // Current value of the fit parameter (S in these slides https://www.snolab.ca/snoplus/private/DocDB/0088/008836/008/AmBe.pdf)
        double scale = params.at("alphanclassibd");

        // Only applies < 3.5 MeV
        if (obs_vals.at(0) < 3.5)
        {
          // Convert energy to bin nuber to get nominal efficiency
          int bin = (obs_vals.at(0) - 1.0) / 0.05;
          double eff = efficiencies[bin];

          // And scale the nominal efficiency accordingly
          double scaledeff = scale * eff;

          // But scaled efficiency must still be between 0 and 1
          if (scaledeff < 0)
            scaledeff = 0;
          else if (scaledeff > 1.0)
            scaledeff = 1.0;

          return scaledeff;
        }
        else
        {
          return 1.0;
        }
      };

      // The bin-by-bin efficiencies scaling for the alpha,n classifier on the proton recoil alpha,n PDF
      ShapeFunction AlphaNClassAlphaN = [](const ParameterDict &params, const std::vector<double> &obs_vals)
      {
        // Nominal efficiency scaling for each bin
        std::vector<double> efficiencies = {0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8,
                                            0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8,
                                            0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8,
                                            0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8,
                                            0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8};

        // Current value of the fit parameter (a in these slides https://www.snolab.ca/snoplus/private/DocDB/0088/008836/008/AmBe.pdf)
        double a = params.at("alphanclassalphan");

        // Only applies < 3.5 MeV
        if (obs_vals.at(0) < 3.5)
        {
          // Convert energy to bin nuber to get nominal efficiency
          int bin = (obs_vals.at(0) - 1.0) / 0.05;
          double eff = efficiencies[bin];

          // And scale the nominal efficiency accordingly
          double scale = 1 + pow(a, 2) * pow(obs_vals.at(0), 3);
          double scaledeff = scale * eff;

          // But scaled efficiency must still be between 0 and 1
          if (scaledeff < 0)
            scaledeff = 0;
          else if (scaledeff > 1.0)
            scaledeff = 1.0;

          return scaledeff;
        }
        else
        {
          return 1.0;
        }
      };

      // Shape function to just scale pdfs up or down on top of their normalisations
      ShapeFunction reac_norm = [](const ParameterDict &params, const std::vector<double> &obs_vals)
      {
        return params.at("reactor_nubar_norm");
      };

      // Shape function to just scale pdfs up or down on top of their normalisations
      ShapeFunction geou_norm = [](const ParameterDict &params, const std::vector<double> &obs_vals)
      {
        return params.at("geonu_U_norm");
      };

      // Shape function to just scale pdfs up or down on top of their normalisations
      ShapeFunction geoth_norm = [](const ParameterDict &params, const std::vector<double> &obs_vals)
      {
        return params.at("geonu_Th_norm");
      };

      // Shape function to just scale pdfs up or down on top of their normalisations
      ShapeFunction alphanpr_norm = [](const ParameterDict &params, const std::vector<double> &obs_vals)
      {
        return params.at("alphan_PRecoil_norm");
      };

      // Shape function to just scale pdfs up or down on top of their normalisations
      ShapeFunction alphancs_norm = [](const ParameterDict &params, const std::vector<double> &obs_vals)
      {
        return params.at("alphan_CScatter_norm");
      };

      // Shape function to just scale pdfs up or down on top of their normalisations
      ShapeFunction alphanoe_norm = [](const ParameterDict &params, const std::vector<double> &obs_vals)
      {
        return params.at("alphan_OExcited_norm");
      };

      // Shape function to just scale pdfs up or down on top of their normalisations
      ShapeFunction bipolike_norm = [](const ParameterDict &params, const std::vector<double> &obs_vals)
      {
        return params.at("bipolike_norm");
      };

      // Shape function to just scale pdfs up or down on top of their normalisations
      ShapeFunction atmospheric_norm = [](const ParameterDict &params, const std::vector<double> &obs_vals)
      {
        return params.at("atmospheric_norm");
      };

      Shape *shape = new Shape(name);

      if (function == "OscProbGrid")
      {
        shape->SetShapeFunction(OscProbGrid, paramnamevec_);
        shape->RenameParameter(paramnamevec_.at(0), "deltam21");
        shape->RenameParameter(paramnamevec_.at(1), "theta12");
        ParameterDict params({{"deltam21", paramvals_[paramnamevec_.at(0)]}, {"theta12", paramvals_[paramnamevec_.at(1)]}});
        shape->SetParameters(params);
      }
      else if (function == "OscProb")
      {
        shape->SetShapeFunction(OscProb, paramnamevec_);
        shape->RenameParameter(paramnamevec_.at(0), "deltam21");
        shape->RenameParameter(paramnamevec_.at(1), "theta12");
        ParameterDict params({{"deltam21", paramvals_[paramnamevec_.at(0)]}, {"theta12", paramvals_[paramnamevec_.at(1)]}});
        shape->SetParameters(params);
      }
      else if (function == "AlphaNClassIBD")
      {
        shape->SetShapeFunction(AlphaNClassIBD, paramnamevec_);
        shape->RenameParameter(paramnamevec_.at(0), "alphanclassibd");
        ParameterDict params({{"alphanclassibd", paramvals_[paramnamevec_.at(0)]}});
        shape->SetParameters(params);
      }
      else if (function == "AlphaNClassAlphaN")
      {
        shape->SetShapeFunction(AlphaNClassAlphaN, paramnamevec_);
        shape->RenameParameter(paramnamevec_.at(0), "alphanclassalphan");
        ParameterDict params({{"alphanclassalphan", paramvals_[paramnamevec_.at(0)]}});
        shape->SetParameters(params);
      }
      else if (function == "reactor_nubar_norm")
      {
        shape->SetShapeFunction(reac_norm, paramnamevec_);
        shape->RenameParameter(paramnamevec_.at(0), "reactor_nubar_norm");
        ParameterDict params({{"reactor_nubar_norm", paramvals_[paramnamevec_.at(0)]}});
        shape->SetParameters(params);
      }
      else if (function == "geonu_U_norm")
      {
        shape->SetShapeFunction(geou_norm, paramnamevec_);
        shape->RenameParameter(paramnamevec_.at(0), "geonu_U_norm");
        ParameterDict params({{"geonu_U_norm", paramvals_[paramnamevec_.at(0)]}});
        shape->SetParameters(params);
      }
      else if (function == "geonu_Th_norm")
      {
        shape->SetShapeFunction(geoth_norm, paramnamevec_);
        shape->RenameParameter(paramnamevec_.at(0), "geonu_Th_norm");
        ParameterDict params({{"geonu_Th_norm", paramvals_[paramnamevec_.at(0)]}});
        shape->SetParameters(params);
      }
      else if (function == "alphan_PRecoil_norm")
      {
        shape->SetShapeFunction(alphanpr_norm, paramnamevec_);
        shape->RenameParameter(paramnamevec_.at(0), "alphan_PRecoil_norm");
        ParameterDict params({{"alphan_PRecoil_norm", paramvals_[paramnamevec_.at(0)]}});
        shape->SetParameters(params);
      }
      else if (function == "alphan_CScatter_norm")
      {
        shape->SetShapeFunction(alphancs_norm, paramnamevec_);
        shape->RenameParameter(paramnamevec_.at(0), "alphan_CScatter_norm");
        ParameterDict params({{"alphan_CScatter_norm", paramvals_[paramnamevec_.at(0)]}});
        shape->SetParameters(params);
      }
      else if (function == "alphan_OExcited_norm")
      {
        shape->SetShapeFunction(alphanoe_norm, paramnamevec_);
        shape->RenameParameter(paramnamevec_.at(0), "alphan_OExcited_norm");
        ParameterDict params({{"alphan_OExcited_norm", paramvals_[paramnamevec_.at(0)]}});
        shape->SetParameters(params);
      }
      else if (function == "bipolike_norm")
      {
        shape->SetShapeFunction(bipolike_norm, paramnamevec_);
        shape->RenameParameter(paramnamevec_.at(0), "bipolike_norm");
        ParameterDict params({{"bipolike_norm", paramvals_[paramnamevec_.at(0)]}});
        shape->SetParameters(params);
      }
      else if (function == "atmospheric_norm")
      {
        shape->SetShapeFunction(atmospheric_norm, paramnamevec_);
        shape->RenameParameter(paramnamevec_.at(0), "atmospheric_norm");
        ParameterDict params({{"atmospheric_norm", paramvals_[paramnamevec_.at(0)]}});
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
