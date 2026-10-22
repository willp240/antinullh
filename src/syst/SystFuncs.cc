#include "SystFuncs.hh"
#include <OscParams.hh>
#include <Exceptions.h>
#include <cmath>
#include <utility>
#include <map>
#include <stdexcept>

namespace antinufit
{

  namespace
  {
    struct OscProbTerms
    {
      double nuEnergy;
      double dEigen10, dEigen20, dEigen21;
      double X10, X20, X21;
      double Evaluate(double baseline) const;
    };

    // The part of OscProb that depends on neutrino energy and the oscillation parameters but not on the
    // baseline. OscProbTerms::Evaluate is numerically equivalent to OscProb, up to floating-point
    // rounding, while allowing the energy-dependent terms to be reused across reactor baselines
    OscProbTerms OscProbEnergyTerms(double nuEnergy, double fDmSqr21, double fSSqrTheta12)
    {
      Double_t nuE_parent = nuEnergy;

      Double_t fDmSqr32 = 2.451e-3;
      Double_t fSSqrTheta13 = 0.0216;
      Double_t Ne = 8.13e23;

      Double_t fDmSqr31 = fDmSqr32 + fDmSqr21;

      Double_t H_ee_vac = fDmSqr21 * (fSSqrTheta12 * (1 - fSSqrTheta13) - (1.0 / 3.0)) + fDmSqr31 * (fSSqrTheta13 - (1.0 / 3.0));
      Double_t H_neq2 = (1 - fSSqrTheta13) * (fDmSqr21 * fDmSqr21 * fSSqrTheta12 * (1 + fSSqrTheta12 * (fSSqrTheta13 - 1)) + fDmSqr31 * fDmSqr31 * fSSqrTheta13 - 2.0 * fDmSqr21 * fDmSqr31 * fSSqrTheta12 * fSSqrTheta13);

      Double_t a0_vac = -(2.0 / 27.0) * (fDmSqr21 * fDmSqr21 * fDmSqr21 + fDmSqr31 * fDmSqr31 * fDmSqr31) + (1.0 / 9.0) * (fDmSqr21 * fDmSqr21 * fDmSqr31 + fDmSqr21 * fDmSqr31 * fDmSqr31);
      Double_t a1_vac = (1.0 / 3.0) * (fDmSqr21 * fDmSqr31 - fDmSqr21 * fDmSqr21 - fDmSqr31 * fDmSqr31);
      Double_t Y_ee_vac = (2.0 / 3.0) * a1_vac + H_ee_vac * H_ee_vac + H_neq2;

      Double_t alpha = -2.535e-31 * Ne; // conversion factor in eV2/MeV
      Double_t A_CC = alpha * nuE_parent; // for A_CC in [eV^2] and nuE in [MeV]
      Double_t alpha_1 = H_ee_vac * A_CC + (1.0 / 3.0) * A_CC * A_CC;

      Double_t a0 = a0_vac - Y_ee_vac * A_CC - (1.0 / 3.0) * H_ee_vac * A_CC * A_CC - (2.0 / 27.0) * A_CC * A_CC * A_CC;
      Double_t a1 = a1_vac - alpha_1;
      Double_t Y_ee = Y_ee_vac + (2.0 / 3.0) * alpha_1;
      Double_t H_ee = H_ee_vac + (2.0 / 3.0) * A_CC;

      // Get eigenvalues of H, and constants X and theta
      Double_t arcCos = (1.0 / 3.0) * acos(1.5 * (a0 / a1) * sqrt(-3.0 / a1));
      Double_t preFact = 2.0 * sqrt(-a1 / 3.0);

      Double_t eigen[3];
      Double_t X[3];
      for (int i = 0; i < 3; ++i)
      {
        eigen[i] = preFact * cos(arcCos - (2.0 / 3.0) * M_PI * i);
        X[i] = (1.0 / 3.0) + (eigen[i] * H_ee + Y_ee) / (3.0 * eigen[i] * eigen[i] + a1);
      }

      OscProbTerms terms;
      terms.nuEnergy = nuE_parent;
      terms.dEigen10 = eigen[1] - eigen[0];
      terms.dEigen20 = eigen[2] - eigen[0];
      terms.dEigen21 = eigen[2] - eigen[1];
      terms.X10 = X[1] * X[0];
      terms.X20 = X[2] * X[0];
      terms.X21 = X[2] * X[1];
      return terms;
    }

    double OscProbTerms::Evaluate(double baseline) const
    {
      Double_t scale = 1.267e3 * baseline / nuEnergy; // for nuE in [MeV] and baseline in [km]
      Double_t s_10 = sin(scale * dEigen10);
      Double_t s_20 = sin(scale * dEigen20);
      Double_t s_21 = sin(scale * dEigen21);
      return 4.0 * (X10 * s_10 * s_10 + X20 * s_20 * s_20 + X21 * s_21 * s_21);
    }

  }

  double BirksLaw(double energy, double birksConstant)
  {
    const double birks_const = 0.074;
    return ((1 + (birks_const * energy)) / (1 + (birksConstant * energy))) * energy;
  }

  double OscProb(double baseline, double nuEnergy, double fDmSqr21, double fSSqrTheta12)
  {
    return OscProbEnergyTerms(nuEnergy, fDmSqr21, fSSqrTheta12).Evaluate(baseline);
  }

  ShapeFunction MakeReactorOscillationShape(const BinnedED &reactorMC,
                                           const std::unordered_map<int, double> &baselines,
                                           const std::string &theta12Name)
  {
    if (reactorMC.GetNDims() != 3 || reactorMC.GetObservables().at(1) != "reactorIndex" ||
        reactorMC.GetObservables().at(2) != "nu_energy")
      throw ValueError("Reactor oscillation needs MC axes (energy, reactorIndex, nu_energy)");
    Theta12ToSinSq(theta12Name, 0.0); // Check the configured angle convention
    const BinAxis energyAxis = reactorMC.GetAxes().GetAxis(0);
    const BinAxis &reactorAxis = reactorMC.GetAxes().GetAxis(1);
    const BinAxis &trueAxis = reactorMC.GetAxes().GetAxis(2);
    const size_t nR = reactorAxis.GetNBins(), nE = trueAxis.GetNBins();

    // Keep only occupied MC bins. Each row belongs to one reconstructed-energy bin;
    // its entries contain a reactor/true-energy index and the MC contribution
    using Entry = std::pair<size_t, double>;
    std::vector<std::vector<Entry>> rows(energyAxis.GetNBins());
    std::vector<bool> occupied(nR * nE, false);
    double total = 0;
    const std::vector<double> counts = reactorMC.GetBinContents();
    for (size_t b = 0; b < counts.size(); ++b)
      if (counts[b] != 0)
      {
        if (!std::isfinite(counts[b]) || counts[b] < 0)
          throw ValueError("Reactor MC must contain finite, nonnegative counts");
        rows[b / (nR * nE)].emplace_back(b % (nR * nE), counts[b]);
        occupied[b % (nR * nE)] = true;
        total += counts[b];
      }
    if (!(total > 0))
      throw ValueError("Reactor oscillation needs nonempty MC");
    std::vector<double> rowSums(rows.size(), 0);
    for (size_t i = 0; i < rows.size(); ++i)
      for (Entry &entry : rows[i])
      {
        entry.second /= total;
        rowSums[i] += entry.second;
      }

    // Group occupied reactor/true-energy pairs by true energy, so the expensive
    // energy-dependent probability terms are calculated once per true-energy bin
    std::vector<std::vector<Entry>> reactors(nE);
    for (size_t r = 0; r < nR; ++r)
    {
      const double low = reactorAxis.GetBinLowEdge(r);
      if (std::fabs(reactorAxis.GetBinHighEdge(r) - low - 1) > 1e-9 ||
          std::fabs(low - std::round(low)) > 1e-9)
        throw ValueError("Reactor-axis bins must be one unit wide and start on an integer");
      for (size_t j = 0; j < nE; ++j)
        if (occupied[r * nE + j])
        {
          const std::unordered_map<int, double>::const_iterator baseline = baselines.find(int(std::lround(low)));
          if (baseline == baselines.end())
            throw ValueError("No baseline for reactor index " + std::to_string(int(std::lround(low))));
          reactors[j].emplace_back(r * nE + j, baseline->second);
        }
    }
    std::vector<double> energies(nE);
    for (size_t j = 0; j < nE; ++j)
    {
      energies[j] = trueAxis.GetBinCentre(j);
    }

    return [energyAxis, rows = std::move(rows), rowSums = std::move(rowSums),
            reactors = std::move(reactors), energies = std::move(energies),
            theta12Name, probabilities = std::vector<double>(nR * nE),
            cached = false, lastDm21 = 0.0, lastTheta12 = 0.0]
           (const ParameterDict &params, const std::vector<double> &binCentre) mutable
    {
      const double dm21 = params.at("deltam21"), theta12 = params.at(theta12Name);
      if (!cached || dm21 != lastDm21 || theta12 != lastTheta12)
      {
        const double sinSq = Theta12ToSinSq(theta12Name, theta12);
        if (!std::isfinite(dm21) || !std::isfinite(theta12) || !std::isfinite(sinSq) || sinSq < 0 || sinSq > 1)
          throw ParameterError("Oscillation parameters must be finite and sin^2(theta12) must be in [0, 1]");
        for (size_t j = 0; j < energies.size(); ++j)
        {
          if (reactors[j].empty())
            continue;
          const OscProbTerms terms = OscProbEnergyTerms(energies[j], dm21, sinSq);
          for (const Entry &reactor : reactors[j])
          {
            probabilities[reactor.first] = 1 - terms.Evaluate(reactor.second);
          }
        }
        lastDm21 = dm21;
        lastTheta12 = theta12;
        cached = true;
      }
      const size_t i = energyAxis.FindBin(binCentre.at(0));
      double weighted = 0;
      for (const Entry &entry : rows[i])
        weighted += entry.second * probabilities[entry.first];
      return rowSums[i] > 0 ? weighted / rowSums[i] : 1.0;
    };
  }

  ShapeFunction MakeAlphaNClassifierShape(const std::string &functionName,
                                          const std::string &parameterName)
  {
    // Nominal efficiencies in 0.05 MeV bins from the classifier inputs
    static const std::map<std::string, std::vector<double>> efficiencies = {
      {"AlphaNClassReacPPO", {0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.860731, 0.861338, 0.869968, 0.873814, 0.873137, 0.873502, 0.879028, 0.883524, 0.884052, 0.882041, 0.886878, 0.889255, 0.891738, 0.896534, 0.898181, 0.903118, 0.906836, 0.905997, 0.909549, 0.912651, 0.913376, 0.913353, 0.919331, 0.918609, 0.921336, 0.923, 0.925586, 0.925844, 0.927308, 0.92664, 0.928934, 0.931587, 0.934654, 0.937111, 0.936353, 0.935488, 0.937217, 0.936887, 0.939664, 0.937311, 0.941545, 0.942837, 0.943906, 0.944927, 0.945262, 0.944521, 0.947274, 0.948066, 0.95056, 0.94877, 0.95, 0.951169, 0.951918, 0.952573, 0.954404, 0.957312, 0.954584, 0.956432, 0.95394, 0.959006, 0.956611, 0.958433, 0.958806, 0.961096, 0.959182, 0.961301, 0.963674, 0.960559, 0.963169, 0.962637, 0.962192, 0.964002, 0.965879, 0.966453, 0.967481, 0.968983, 0.969748, 0.967684, 0.969748, 0.968921, 0.969408, 0.96981, 0.967608, 0.968713, 0.97169, 0.970209, 0.972268, 0.974546, 0.970594, 0.974388, 0.971034, 0.973232, 0.972703, 0.973704, 0.974172, 0.974487, 0.971795, 0.9747, 0.975558, 0.97901, 0.97926, 0.979101, 0.979386, 0.980533, 0.975214, 0.978767, 0.976117, 0.977475, 0.978156, 0.979894, 0.981468, 0.980643, 0.983536, 0.984319, 0.984496, 0.982833, 0.979479, 0.982598, 0.985887, 0.986774, 0.986689, 0.983193, 0.984095, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}},
      {"AlphaNClassReacBisMSB", {0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.806548, 0.835778, 0.860114, 0.860511, 0.861445, 0.866607, 0.872169, 0.875728, 0.879387, 0.880189, 0.885962, 0.884388, 0.890787, 0.896026, 0.896036, 0.900734, 0.903249, 0.90666, 0.908125, 0.911743, 0.912975, 0.916396, 0.919068, 0.921901, 0.924725, 0.927896, 0.928439, 0.928985, 0.932167, 0.934756, 0.936884, 0.936454, 0.939462, 0.942235, 0.943137, 0.946213, 0.947471, 0.950546, 0.95214, 0.952979, 0.955625, 0.954861, 0.958558, 0.960074, 0.959584, 0.961882, 0.963916, 0.963997, 0.965395, 0.965859, 0.968, 0.969706, 0.970332, 0.972599, 0.972303, 0.973358, 0.973549, 0.975858, 0.974792, 0.976407, 0.978204, 0.9787, 0.978272, 0.980631, 0.98141, 0.981719, 0.983176, 0.983758, 0.982372, 0.9837, 0.985131, 0.985506, 0.987578, 0.986775, 0.98647, 0.986375, 0.987629, 0.988508, 0.989703, 0.988062, 0.989839, 0.990823, 0.990985, 0.9906, 0.990882, 0.990478, 0.990748, 0.992284, 0.992133, 0.993124, 0.992274, 0.993092, 0.992593, 0.992843, 0.99282, 0.99424, 0.994556, 0.993289, 0.994567, 0.994571, 0.995203, 0.994617, 0.995686, 0.994558, 0.994818, 0.995894, 0.994399, 0.996132, 0.996935, 0.996536, 0.994534, 0.997692, 0.996349, 0.996895, 0.996946, 0.995282, 0.996984, 0.997012, 0.996868, 0.99804, 0.998264, 0.997658, 0.999455, 0.998093, 0.999265, 0.998221, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}},
      {"AlphaNClassGeoUPPO", {0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.847601, 0.862521, 0.869404, 0.872788, 0.868274, 0.874179, 0.875593, 0.882563, 0.880354, 0.891766, 0.888742, 0.887482, 0.891711, 0.892535, 0.899591, 0.897328, 0.904896, 0.901258, 0.900262, 0.91033, 0.911279, 0.912267, 0.916414, 0.92028, 0.924047, 0.918843, 0.920548, 0.925015, 0.928582, 0.933845, 0.929945, 0.926582, 0.931034, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}},
      {"AlphaNClassGeoUBisMSB", {0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.77118, 0.841852, 0.861376, 0.865902, 0.861295, 0.870079, 0.87114, 0.876489, 0.881636, 0.884859, 0.88409, 0.883534, 0.880806, 0.880433, 0.895535, 0.897326, 0.9011, 0.905386, 0.90821, 0.916148, 0.914288, 0.917934, 0.914414, 0.918733, 0.919882, 0.927415, 0.930227, 0.932457, 0.932278, 0.932484, 0.934144, 0.935567, 0.925071, 0.920193, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}},
      {"AlphaNClassGeoThPPO", {0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.869974, 0.8654, 0.869456, 0.871193, 0.872548, 0.874571, 0.868029, 0.880168, 0.882073, 0.881392, 0.880617, 0.886276, 0.879296, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}},
      {"AlphaNClassGeoThBisMSB", {0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.842345, 0.859232, 0.863687, 0.867973, 0.873233, 0.880205, 0.883338, 0.878529, 0.885726, 0.888816, 0.875962, 0.874857, 0.845707, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}},
      {"AlphaNClassBPLikePPO", {0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.847718, 0.853238, 0.857366, 0.864198, 0.86795, 0.872599, 0.875782, 0.879574, 0.881412, 0.887901, 0.887572, 0.891668, 0.89631, 0.89612, 0.900281, 0.901399, 0.906291, 0.906306, 0.905812, 0.908573, 0.909191, 0.911198, 0.912269, 0.914091, 0.915853, 0.916326, 0.917088, 0.919463, 0.919963, 0.921576, 0.922873, 0.923135, 0.925275, 0.927401, 0.927655, 0.929038, 0.930798, 0.932076, 0.932964, 0.933235, 0.933717, 0.935461, 0.937579, 0.937692, 0.938867, 0.939073, 0.939769, 0.939569, 0.940584, 0.941851, 0.941465, 0.94102, 0.939823, 0.941513, 0.936784, 0.94052, 0.929245, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}},
      {"AlphaNClassBPLikeBisMSB", {0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.825317, 0.831639, 0.836172, 0.839495, 0.842258, 0.84561, 0.851943, 0.857137, 0.859025, 0.86411, 0.867047, 0.870545, 0.873277, 0.876216, 0.881272, 0.88253, 0.886845, 0.889446, 0.890827, 0.893902, 0.895305, 0.898821, 0.903277, 0.908194, 0.912053, 0.914033, 0.916743, 0.919293, 0.921249, 0.923123, 0.925891, 0.928224, 0.93094, 0.932513, 0.934881, 0.937013, 0.938612, 0.940964, 0.943078, 0.945012, 0.946747, 0.948388, 0.949747, 0.951458, 0.953077, 0.953776, 0.954598, 0.955318, 0.957422, 0.957053, 0.954817, 0.953522, 0.950674, 0.944864, 0.939431, 0.934121, 0.932976, 0.919776, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}},
      {"AlphaNClassAlphaNPPO", {0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.379924, 0.36871, 0.346533, 0.334301, 0.321632, 0.307839, 0.29835, 0.282151, 0.275494, 0.263849, 0.251177, 0.243713, 0.237546, 0.230593, 0.218946, 0.213989, 0.206642, 0.20184, 0.189907, 0.186113, 0.181831, 0.172919, 0.167854, 0.16315, 0.157402, 0.152269, 0.151937, 0.147145, 0.143054, 0.137624, 0.135609, 0.13182, 0.126017, 0.125131, 0.123167, 0.121317, 0.119797, 0.11557, 0.116512, 0.114151, 0.110107, 0.110443, 0.112288, 0.108253, 0.107296, 0.108336, 0.111337, 0.110101, 0.110948, 0.107593, 0.102833, 0.105871, 0.104421, 0.107876, 0.111511, 0.110437, 0.0986493, 0.0979467, 0.109627, 0.110945, 0.110842, 0.110141, 0.12345, 0.123803, 0.114362, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.932749, 0.934783, 0.940553, 0.933628, 0.934715, 0.939334, 0.927531, 0.929472, 0.920963, 0.917469, 0.924277, 0.914058, 0.909152, 0.905714, 0.897545, 0.898317, 0.900198, 1, 1, 1, 1, 1, 1, 0.976, 0.982552, 0.981351, 0.983996, 0.984207, 0.983028, 0.983027, 0.981115, 0.977986, 0.976765, 0.97374, 0.972434, 0.972522, 0.971942, 0.966004, 0.96786, 0.966127, 0.964405, 0.970949, 0.962814, 0.958261, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}},
      {"AlphaNClassAlphaNBisMSB", {0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.22799, 0.211668, 0.203507, 0.195213, 0.187173, 0.172563, 0.163788, 0.154083, 0.148568, 0.14172, 0.12975, 0.124008, 0.118391, 0.111164, 0.107499, 0.100328, 0.0966493, 0.0912828, 0.0855004, 0.0833023, 0.0792089, 0.0763367, 0.0717646, 0.0694904, 0.068143, 0.0637377, 0.0619927, 0.0617361, 0.0595614, 0.0562207, 0.0552019, 0.0548886, 0.0542087, 0.0523142, 0.0512037, 0.0519437, 0.0500508, 0.0493043, 0.049204, 0.0491005, 0.0480187, 0.0494916, 0.0483119, 0.044915, 0.0466631, 0.0487034, 0.0484837, 0.0467251, 0.0458571, 0.0453976, 0.0478289, 0.049885, 0.0501955, 0.0531463, 0.051851, 0.0498724, 0.049746, 0.0488046, 0.0553527, 0.0557996, 0.0606216, 0.0675039, 0.0635338, 0.0597531, 0.0832274, 0.0797774, 1, 1, 1, 1, 1, 1, 1, 0.952102, 0.965794, 0.972275, 0.978573, 0.981465, 0.976859, 0.976728, 0.97601, 0.977394, 0.979429, 0.973397, 0.97174, 0.969115, 0.969677, 0.969263, 0.96476, 0.956171, 0.968675, 1, 1, 1, 1, 1, 1, 0.973254, 0.977254, 0.986588, 0.992719, 0.993739, 0.994637, 0.995701, 0.996282, 0.996361, 0.995705, 0.996196, 0.995251, 0.995545, 0.994787, 0.99464, 0.993001, 0.992579, 0.993831, 0.992841, 0.98965, 0.986787, 0.992063, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}},
    };
    const std::map<std::string, std::vector<double>>::const_iterator table =
        efficiencies.find(functionName);
    if (table == efficiencies.end())
      throw std::invalid_argument("Unknown alpha,n classifier shape: " + functionName);
    const std::vector<double> *values = &table->second;
    const bool alphaN = functionName == "AlphaNClassAlphaNPPO" ||
                        functionName == "AlphaNClassAlphaNBisMSB";
    return [values, parameterName, alphaN](const ParameterDict &params,
                                            const std::vector<double> &obsVals)
    {
      const double energy = obsVals.at(0);
      if (energy >= 3.5)
        return 1.0;
      const size_t bin = static_cast<size_t>(energy / 0.05);
      const double parameter = params.at(parameterName);
      const double scale = alphaN ? 1.0 + parameter * parameter * energy * energy * energy
                                  : parameter;
      return std::max(0.0, std::min(1.0, scale * values->at(bin)));
    };
  }

}
