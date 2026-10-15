#pragma once

// OXO headers
#include <Systematic.h>
#include <BinnedED.h>
#include <DataSet.h>

// c++ headers
#include <string>
#include <vector>
#include <unordered_map>
#include <cstddef>
#include <set>

namespace antinufit
{
  /*
   * Reactor oscillation applied to a 1D reconstructed-energy PDF.
   *
   * The unoscillated reactor MC is kept as a sparse 3D histogram H[i][r][j]
   * (reconstructed energy i, reactor index r, true antineutrino energy j). Each Construct() computes
   *
   *   P[r][j] = 1 - OscProb(L_r, E_j, deltam21, sin^2(theta12))
   *   R[i]    = sum_{r,j} H[i][r][j] / N_R * P[r][j]
   *
   * N_R is the total reactor MC count in the source histogram, not the number of reactor cores.
   * BinCentre evaluates each true-energy bin at its centre, approximating event-level weighting.
   * The response is diagonal: w_i = R[i] / H1D[i], where H1D[i] = sum_{r,j} H[i][r][j] / N_R is the
   * registered 1D PDF. Applied to that PDF it returns exactly R, so the reactor PDF integral becomes
   * the MC-weighted mean survival probability.
   *
   * The PDF this acts on must be built from the same events as the 3D source, and no systematic that
   * moves content between energy bins may be applied before it.
   */
  class OscWeight1D : public Systematic
  {
  public:
    enum ProbMode
    {
      BinCentre, // P evaluated at the true-energy bin centre
      BinAverage // P averaged over the true-energy bin, assuming the phase is linear in E across the bin
    };

    // Angle name selects units: theta12 (degrees), sintheta12, or sinsqtheta12.
    OscWeight1D(const std::string &name_,
                const std::string &dm21Name_,
                const std::string &theta12Name_,
                const std::unordered_map<int, double> &indexDistance_,
                ProbMode mode_ = BinCentre);

    // Setup, once per reactor PDF: h3d_ and h1d_ are unnormalised MC counts, before padding.
    // Throws unless every row of h3d_ sums to the matching bin of h1d_.
    void SetSource(const BinnedED &h3d_, const BinnedED &h1d_);
    bool HasSource() const { return fNR > 0; }

    // Setup check on the MC events: every reactorIndex must be an integer inside the reactor axis and
    // have a baseline. Returns the number of events whose true energy is outside the true-energy axis
    // (they are clamped into its edge bins).
    static size_t CheckEvents(DataSet *data_, const AxisCollection &axes3d_,
                              const std::unordered_map<int, double> &indexDistance_);

    void Construct();

    // Summaries of the source and of the last Construct()
    double GetNR() const { return fNR; }
    size_t GetNNonZero() const { return fVal.size(); }
    size_t GetNPairs() const { return fPairL.size(); }
    size_t GetN3DBins() const { return fN3DBins; }
    const std::vector<double> &GetR() const { return fR; }
    const std::vector<double> &GetSourceMarginal() const { return fRowSum; }
    double GetMeanSurvival() const;
    ProbMode GetMode() const { return fMode; }

    // Accumulated wall time in Construct()
    double GetConstructSeconds() const { return fConstructSeconds; }
    size_t GetNConstruct() const { return fNConstruct; }

    // FitComponent interface
    void SetParameter(const std::string &name_, double value_);
    double GetParameter(const std::string &name_) const;
    void SetParameters(const ParameterDict &params_);
    ParameterDict GetParameters() const;
    size_t GetParameterCount() const { return 2; }
    std::set<std::string> GetParameterNames() const;
    void RenameParameter(const std::string &old_, const std::string &new_);
    std::string GetName() const { return fName; }
    void SetName(const std::string &name_) { fName = name_; }

    static const std::string kReactorBranch;
    static const std::string kNuEnergyBranch;

  private:
    std::string fName;
    std::string fDm21Name;
    std::string fTheta12Name;
    double fDm21;
    double fTheta12Value; // Value in the configured representation; converted only for evaluation.
    std::string fTheta12Convention; // Preserved when the exposed parameter is renamed.
    std::unordered_map<int, double> fIndexDistance;
    ProbMode fMode;

    // Sparse source, grouped by reconstructed-energy bin i: entries fRowStart[i] .. fRowStart[i+1]-1
    std::vector<size_t> fRowStart;
    std::vector<double> fVal;       // H[i][r][j] / N_R
    std::vector<unsigned> fPairIdx; // index of (r, j) in the occupied-pair list
    std::vector<double> fRowSum;    // sum_{r,j} H[i][r][j] / N_R, the normalised H1D

    // Occupied (r, j) pairs, grouped by true-energy bin j: pairs fJStart[j] .. fJStart[j+1]-1
    std::vector<size_t> fJStart;
    std::vector<double> fPairL;       // baseline of each pair's reactor [km]
    std::vector<double> fEnergyCentre; // true-energy bin centres [MeV]
    std::vector<double> fEnergyWidth;  // true-energy bin widths [MeV]

    std::vector<double> fP; // P_ee per occupied pair
    std::vector<double> fR;
    std::vector<unsigned long long> fDiagIdx;

    double fNR;
    size_t fN3DBins;
    double fConstructSeconds;
    size_t fNConstruct;
  };
}
