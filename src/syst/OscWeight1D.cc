#include <OscWeight1D.hh>
#include <Functions.hh>
#include <OscParams.hh>

// OXO headers
#include <Event.h>
#include <Exceptions.h>

// c++ headers
#include <chrono>
#include <algorithm>
#include <cmath>
#include <sstream>
#include <stdexcept>

namespace antinufit
{
  const std::string OscWeight1D::kReactorBranch = "reactorIndex";
  const std::string OscWeight1D::kNuEnergyBranch = "nu_energy";

  namespace
  {
    // Reactor index of a reactor-axis bin. Bins must be one unit wide and start on an integer, so that
    // each bin holds exactly one reactor index.
    int ReactorIndexOfBin(const BinAxis &axis_, size_t bin_)
    {
      const double low = axis_.GetBinLowEdge(bin_);
      const double high = axis_.GetBinHighEdge(bin_);
      if (std::fabs(high - low - 1) > 1e-9 || std::fabs(low - std::round(low)) > 1e-9)
      {
        std::ostringstream msg;
        msg << "OscWeight1D: reactor axis bin " << bin_ << " is [" << low << ", " << high
            << "); every bin must be one unit wide and start on an integer";
        throw std::runtime_error(msg.str());
      }
      return int(std::lround(low));
    }

    // <sin^2(phi)> over a true-energy bin, for phi linear in E across the bin with total spread w
    inline double MeanSinSq(double phi_, double w_)
    {
      const double s = sin(phi_);
      const double sinc = (w_ < 1e-4) ? 1 - w_ * w_ / 6 : sin(w_) / w_;
      return s * s * sinc + 0.5 * (1 - sinc);
    }
  }

  OscWeight1D::OscWeight1D(const std::string &name_,
                           const std::string &dm21Name_,
                           const std::string &theta12Name_,
                           const std::unordered_map<int, double> &indexDistance_,
                           ProbMode mode_)
      : fName(name_), fDm21Name(dm21Name_), fTheta12Name(theta12Name_), fDm21(0), fTheta12Value(0), fTheta12Convention(theta12Name_),
        fIndexDistance(indexDistance_), fMode(mode_), fNR(0), fN3DBins(0),
        fConstructSeconds(0), fNConstruct(0)
  {
    // Validate the representation before any MC is loaded.
    Theta12ToSinSq(fTheta12Convention, 0.0);
  }

  size_t OscWeight1D::CheckEvents(DataSet *data_, const AxisCollection &axes3d_,
                                  const std::unordered_map<int, double> &indexDistance_)
  {
    const BinAxis &reactorAxis = axes3d_.GetAxis(1);
    const BinAxis &nuEAxis = axes3d_.GetAxis(2);
    size_t nClamped = 0;
    for (size_t iEv = 0; iEv < data_->GetNEntries(); iEv++)
    {
      Event ev = data_->GetEntry(iEv);
      const double reactor = ev.GetDatum(kReactorBranch);
      const double nuE = ev.GetDatum(kNuEnergyBranch);
      const long idx = std::lround(reactor);
      if (std::fabs(reactor - idx) > 1e-6 || reactor < reactorAxis.GetMin() || reactor >= reactorAxis.GetMax())
      {
        std::ostringstream msg;
        msg << "OscWeight1D: event " << iEv << " has " << kReactorBranch << " = " << reactor
            << ", not an integer inside the reactor axis [" << reactorAxis.GetMin() << ", " << reactorAxis.GetMax() << ")";
        throw std::runtime_error(msg.str());
      }
      if (indexDistance_.find(int(idx)) == indexDistance_.end())
      {
        std::ostringstream msg;
        msg << "OscWeight1D: event " << iEv << " has " << kReactorBranch << " = " << idx
            << ", which has no baseline in the reactor JSON";
        throw std::runtime_error(msg.str());
      }
      if (nuE < nuEAxis.GetMin() || nuE >= nuEAxis.GetMax())
        nClamped++;
    }
    return nClamped;
  }

  void OscWeight1D::SetSource(const BinnedED &h3d_, const BinnedED &h1d_)
  {
    if (HasSource())
      throw std::runtime_error("OscWeight1D " + fName + ": source already set; one instance serves one reactor PDF");
    if (h3d_.GetNDims() != 3 || h1d_.GetNDims() != 1)
      throw std::runtime_error("OscWeight1D " + fName + ": need a 3D source and a 1D PDF");

    const std::vector<std::string> &obs = h3d_.GetObservables();
    if (obs.at(1) != kReactorBranch || obs.at(2) != kNuEnergyBranch)
      throw std::runtime_error("OscWeight1D " + fName + ": 3D source axes must be (energy, " + kReactorBranch + ", " + kNuEnergyBranch + ")");

    const AxisCollection &axes = h3d_.GetAxes();
    const BinAxis &energyAxis = axes.GetAxis(0);
    const BinAxis &reactorAxis = axes.GetAxis(1);
    const BinAxis &nuEAxis = axes.GetAxis(2);
    if (!(energyAxis == h1d_.GetAxes().GetAxis(0)))
      throw std::runtime_error("OscWeight1D " + fName + ": 3D source and 1D PDF have different energy binning");
    if (fAxes.GetNBins() != energyAxis.GetNBins())
      throw std::runtime_error("OscWeight1D " + fName + ": SetAxes() must be given the 1D energy axis of the PDF");

    const size_t nI = energyAxis.GetNBins();
    const size_t nR = reactorAxis.GetNBins();
    const size_t nJ = nuEAxis.GetNBins();
    fN3DBins = nI * nR * nJ;

    // Pass 1: which (r, j) pairs are occupied, and the row sums for the check against H1D
    const std::vector<double> contents = h3d_.GetBinContents();
    std::vector<char> pairUsed(nR * nJ, 0);
    std::vector<double> rowCounts(nI, 0);
    double total = 0;
    size_t nnz = 0;
    for (size_t b = 0; b < fN3DBins; b++)
    {
      const double c = contents[b];
      if (c == 0)
        continue;
      rowCounts[b / (nR * nJ)] += c;
      pairUsed[b % (nR * nJ)] = 1;
      total += c;
      nnz++;
    }

    for (size_t i = 0; i < nI; i++)
    {
      const double h1d = h1d_.GetBinContent(i);
      if (std::fabs(rowCounts[i] - h1d) > 1e-9 * std::max(1.0, h1d))
      {
        std::ostringstream msg;
        msg << "OscWeight1D " << fName << ": energy bin " << i << " holds " << rowCounts[i]
            << " events in the 3D source but " << h1d << " in the 1D PDF";
        throw std::runtime_error(msg.str());
      }
    }
    if (total <= 0)
      throw std::runtime_error("OscWeight1D " + fName + ": empty source");

    // Occupied pairs, grouped by true-energy bin, with their baselines
    std::vector<int> pairIndex(nR * nJ, -1);
    fJStart.assign(nJ + 1, 0);
    fPairL.clear();
    for (size_t j = 0; j < nJ; j++)
    {
      fJStart[j] = fPairL.size();
      for (size_t r = 0; r < nR; r++)
      {
        if (!pairUsed[r * nJ + j])
          continue;
        const int reactor = ReactorIndexOfBin(reactorAxis, r);
        std::unordered_map<int, double>::const_iterator it = fIndexDistance.find(reactor);
        if (it == fIndexDistance.end())
        {
          std::ostringstream msg;
          msg << "OscWeight1D " << fName << ": reactor index " << reactor << " is in the MC but has no baseline";
          throw std::runtime_error(msg.str());
        }
        pairIndex[r * nJ + j] = fPairL.size();
        fPairL.push_back(it->second);
      }
    }
    fJStart[nJ] = fPairL.size();

    fEnergyCentre.resize(nJ);
    fEnergyWidth.resize(nJ);
    for (size_t j = 0; j < nJ; j++)
    {
      fEnergyCentre[j] = nuEAxis.GetBinCentre(j);
      fEnergyWidth[j] = nuEAxis.GetBinHighEdge(j) - nuEAxis.GetBinLowEdge(j);
    }

    // Pass 2: the sparse rows, normalised by N_R
    fNR = total;
    fRowStart.assign(nI + 1, 0);
    fVal.clear();
    fPairIdx.clear();
    fVal.reserve(nnz);
    fPairIdx.reserve(nnz);
    fRowSum.assign(nI, 0);
    for (size_t i = 0; i < nI; i++)
    {
      fRowStart[i] = fVal.size();
      for (size_t rj = 0; rj < nR * nJ; rj++)
      {
        const double c = contents[i * nR * nJ + rj];
        if (c == 0)
          continue;
        fVal.push_back(c / fNR);
        fPairIdx.push_back(pairIndex[rj]);
        fRowSum[i] += c / fNR;
      }
    }
    fRowStart[nI] = fVal.size();

    fP.assign(fPairL.size(), 1);
    fR.assign(nI, 0);
    fDiagIdx.resize(nI);
    for (size_t i = 0; i < nI; i++)
      fDiagIdx[i] = i;
  }

  void OscWeight1D::Construct()
  {
    if (!HasSource())
      throw std::runtime_error("OscWeight1D " + fName + ": Construct() called before SetSource()");
    const double sinSqTheta12 = Theta12ToSinSq(fTheta12Convention, fTheta12Value);
    if (!std::isfinite(fDm21) || !std::isfinite(fTheta12Value) ||
        !std::isfinite(sinSqTheta12) || sinSqTheta12 < 0.0 || sinSqTheta12 > 1.0)
      throw ParameterError("OscWeight1D " + fName +
                           ": oscillation parameters must be finite and sin^2(theta12) must be in [0, 1]");
    const std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();

    // P_ee for every occupied (r, j): the energy-dependent terms once per j, then one pass over reactors
    const size_t nJ = fEnergyCentre.size();
    const double *pairL = fPairL.data();
    double *P = fP.data();
    for (size_t j = 0; j < nJ; j++)
    {
      if (fJStart[j] == fJStart[j + 1])
        continue;
      const OscProbTerms terms = OscProbEnergyTerms(fEnergyCentre[j], fDm21, sinSqTheta12);
      if (fMode == BinCentre)
      {
        for (size_t p = fJStart[j]; p < fJStart[j + 1]; p++)
          P[p] = 1.0 - terms.Evaluate(pairL[p]);
      }
      else
      {
        const double relWidth = fEnergyWidth[j] / fEnergyCentre[j];
        for (size_t p = fJStart[j]; p < fJStart[j + 1]; p++)
        {
          const double scale = 1.267e3 * pairL[p] / terms.nuEnergy;
          const double phi10 = scale * terms.dEigen10;
          const double phi20 = scale * terms.dEigen20;
          const double phi21 = scale * terms.dEigen21;
          P[p] = 1.0 - 4.0 * (terms.X10 * MeanSinSq(phi10, std::fabs(phi10) * relWidth) +
                              terms.X20 * MeanSinSq(phi20, std::fabs(phi20) * relWidth) +
                              terms.X21 * MeanSinSq(phi21, std::fabs(phi21) * relWidth));
        }
      }
    }

    // R_i = sum_{r,j} h[i][r][j] P[r][j], and the diagonal weight that turns the registered PDF into R
    const size_t nI = fR.size();
    const double *val = fVal.data();
    const unsigned *pairIdx = fPairIdx.data();
    std::vector<double> weights(nI, 1);
    for (size_t i = 0; i < nI; i++)
    {
      double sum = 0;
      for (size_t k = fRowStart[i]; k < fRowStart[i + 1]; k++)
        sum += val[k] * P[pairIdx[k]];
      fR[i] = sum;
      // An empty row keeps weight 1, so the PDF's padding in that bin is left alone
      if (fRowSum[i] > 0)
        weights[i] = sum / fRowSum[i];
    }
    fResponse.SetComponents(fDiagIdx, fDiagIdx, weights);

    fConstructSeconds += std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
    fNConstruct++;
  }

  double OscWeight1D::GetMeanSurvival() const
  {
    double sum = 0;
    for (size_t i = 0; i < fR.size(); i++)
      sum += fR[i];
    return sum;
  }

  void OscWeight1D::SetParameter(const std::string &name_, double value_)
  {
    if (name_ == fDm21Name)
      fDm21 = value_;
    else if (name_ == fTheta12Name)
      fTheta12Value = value_;
    else
      throw ParameterError("OscWeight1D: can't set " + name_);
  }

  double OscWeight1D::GetParameter(const std::string &name_) const
  {
    if (name_ == fDm21Name)
      return fDm21;
    if (name_ == fTheta12Name)
      return fTheta12Value;
    throw ParameterError("OscWeight1D: parameter " + name_ + " does not exist");
  }

  void OscWeight1D::SetParameters(const ParameterDict &params_)
  {
    ParameterDict::const_iterator dm21 = params_.find(fDm21Name);
    ParameterDict::const_iterator s12sq = params_.find(fTheta12Name);
    if (dm21 != params_.end())
      fDm21 = dm21->second;
    if (s12sq != params_.end())
      fTheta12Value = s12sq->second;
  }

  ParameterDict OscWeight1D::GetParameters() const
  {
    ParameterDict params;
    params[fDm21Name] = fDm21;
    params[fTheta12Name] = fTheta12Value;
    return params;
  }

  std::set<std::string> OscWeight1D::GetParameterNames() const
  {
    std::set<std::string> names;
    names.insert(fDm21Name);
    names.insert(fTheta12Name);
    return names;
  }

  void OscWeight1D::RenameParameter(const std::string &old_, const std::string &new_)
  {
    if (old_ == fDm21Name)
      fDm21Name = new_;
    else if (old_ == fTheta12Name)
      fTheta12Name = new_;
    else
      throw ParameterError("OscWeight1D: can't rename " + old_);
  }
}
