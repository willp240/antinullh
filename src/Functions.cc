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

    Double_t nuE_parent;
    Double_t nuE;

    RAT::DBLinkPtr linkdb;
    RAT::DB *db = RAT::DB::Get();
    RAT::DS::Run run;
    db->SetAirplaneModeStatus(true);
    db->LoadDefaults();

    linkdb = db->GetLink("OSCILLATIONS");
    Double_t fDmSqr21 = params.at("deltam21");
    Double_t fDmSqr32 = linkdb->GetD("sinsqrtheta13");
    Double_t fSSqrTheta12 = sin(params.at("theta12")) * sin(params.at("theta12"));
    Double_t fSSqrTheta13 = linkdb->GetD("sinsqrtheta13");

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

    Double_t baseline = obs_vals.at(1);

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
    // Calc osc prob as a function of params["theta12"], params["delatm21"], obs_val (energy)

    return fOscProb;
  }

  std::map<std::string, FunctionVariant> functionMap = {
      {"BirksLaw", BirksLaw},
      {"OscProb", OscProb}};

}
