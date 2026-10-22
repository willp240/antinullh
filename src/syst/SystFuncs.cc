#include "SystFuncs.hh"
#include <OscParams.hh>
#include <Exceptions.h>
#include <cmath>
#include <utility>

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

  double BirksLaw(const ParameterDict &params, const double &obs_val)
  {
    const double birks_const = 0.074;
    return ((1 + (birks_const * obs_val)) / (1 + (params.at("birks_constant") * obs_val))) * obs_val;
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

}
