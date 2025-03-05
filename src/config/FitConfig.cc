#include "FitConfig.hh"

namespace antinufit
{

  ParameterDict
  FitConfig::GetMinima() const
  {
    return fMinima;
  }

  ParameterDict
  FitConfig::GetMaxima() const
  {
    return fMaxima;
  }

  ParameterDict
  FitConfig::GetNominals() const
  {
    return fNominals;
  }

  ParameterDict
  FitConfig::GetFakeDataVals() const
  {
    return fFakeDataVals;
  }

  ParameterDict
  FitConfig::GetSigmas() const
  {
    return fSigmas;
  }

  ParameterDict
  FitConfig::GetNBins() const
  {
    return fNbins;
  }

  std::map<std::string, std::string>
  FitConfig::GetTexLabels() const
  {
    return fTexLabels;
  }

  int FitConfig::GetIterations() const
  {
    return fIterations;
  }

  void
  FitConfig::SetIterations(int it_)
  {
    fIterations = it_;
  }

  int FitConfig::GetBurnIn() const
  {
    return fBurnIn;
  }

  void
  FitConfig::SetBurnIn(int b_)
  {
    fBurnIn = b_;
  }

  int FitConfig::GetHMCIterations() const
  {
    return fHMCIterations;
  }

  void
  FitConfig::SetHMCIterations(int it_)
  {
    fHMCIterations = it_;
  }

  int FitConfig::GetHMCBurnIn() const
  {
    return fHMCBurnIn;
  }

  void
  FitConfig::SetHMCBurnIn(int b_)
  {
    fHMCBurnIn = b_;
  }

  double
  FitConfig::GetSigmaScale() const
  {
    return fSigmaScale;
  }

  void
  FitConfig::SetSigmaScale(double b_)
  {
    fSigmaScale = b_;
  }

  bool
  FitConfig::GetBeestonBarlow() const
  {
    return fBeestonBarlow;
  }

  void
  FitConfig::SetBeestonBarlow(bool b_)
  {
    fBeestonBarlow = b_;
  }

  bool
  FitConfig::GetAsimov() const
  {
    return fAsimov;
  }

  void
  FitConfig::SetAsimov(bool b_)
  {
    fAsimov = b_;
  }

  bool
  FitConfig::GetFakeData() const
  {
    return fFakeData;
  }

  void
  FitConfig::SetFakeData(bool b_)
  {
    fFakeData = b_;
  }

  std::string
  FitConfig::GetDatafile() const
  {
    return fDatafile;
  }

  void
  FitConfig::SetDatafile(std::string s_)
  {
    fDatafile = s_;
  }

  double
  FitConfig::GetLivetime() const
  {
    return fLivetime;
  }

  void
  FitConfig::SetLivetime(double d_)
  {
    fLivetime = d_;
  }

  std::set<std::string>
  FitConfig::GetParamNames() const
  {
    return ContainerTools::GetKeys(fMinima);
  }

  const std::string &
  FitConfig::GetOutDir() const
  {
    return fOutDir;
  }

  void
  FitConfig::SetOutDir(const std::string &s_)
  {
    fOutDir = s_;
  }

  int FitConfig::GetNSteps() const
  {
    return fNsteps;
  }

  void
  FitConfig::SetNSteps(int n_)
  {
    fNsteps = n_;
  }

  double
  FitConfig::GetEpsilon() const
  {
    return fEpsilon;
  }

  void
  FitConfig::SetEpsilon(double e_)
  {
    fEpsilon = e_;
  }

  ParameterDict
  FitConfig::GetConstrMeans() const
  {
    return fConstrMeans;
  }

  ParameterDict
  FitConfig::GetConstrSigmas() const
  {
    return fConstrSigmas;
  }

  ParameterDict
  FitConfig::GetConstrRatioMeans() const
  {
    return fConstrRatioMeans;
  }

  ParameterDict
  FitConfig::GetConstrRatioSigmas() const
  {
    return fConstrRatioSigmas;
  }

  std::map<std::string, std::string>
  FitConfig::GetConstrRatioParName() const
  {
    return fConstrRatioParName;
  }

  ParameterDict
  FitConfig::GetConstrCorrs() const
  {
    return fConstrCorrs;
  }

  std::map<std::string, std::string>
  FitConfig::GetConstrCorrParName() const
  {
    return fConstrCorrParName;
  }

  bool
  FitConfig::GetSaveOutputs() const
  {
    return fSaveOutputs;
  }

  void
  FitConfig::SetSaveOutputs(bool b_)
  {
    fSaveOutputs = b_;
  }

  void
  FitConfig::AddParameter(const std::string &name_, double nom_, double min_, double max_, double sigma_, int nbins_, double fdvalue_,
                          std::string label_, double constrMean_, double constrSigma_)
  {

    fConstrMeans[name_] = constrMean_;
    fConstrSigmas[name_] = constrSigma_;

    AddParameter(name_, nom_, min_, max_, sigma_, nbins_, fdvalue_, label_);
  }

  void
  FitConfig::AddParameter(const std::string &name_, double nom_, double min_, double max_, double sigma_, int nbins_, double fdvalue_,
                          std::string label_)
  {
    fMinima[name_] = min_;
    fMaxima[name_] = max_;
    fNominals[name_] = nom_;
    fFakeDataVals[name_] = fdvalue_;
    fSigmas[name_] = sigma_;
    fNbins[name_] = nbins_;
    fTexLabels[name_] = label_;
  }

  void
  FitConfig::AddParameter(const std::string &name_, double nom_, double min_, double max_, double sigma_, int nbins_, double fdvalue_,
                          std::string label_, double constrRatioMean_, double constrRatioSigma_, std::string constrRatioParName_)
  {

    fConstrRatioMeans[name_] = constrRatioMean_;
    fConstrRatioSigmas[name_] = constrRatioSigma_;
    fConstrRatioParName[name_] = constrRatioParName_;

    AddParameter(name_, nom_, min_, max_, sigma_, nbins_, fdvalue_, label_);
  }

  void
  FitConfig::AddParameter(const std::string &name_, double nom_, double min_, double max_, double sigma_, int nbins_, double fdvalue_,
                          std::string label_, double constrMean_, double constrSigma_, std::string constrCorrParName_, double constrCorr_)
  {

    fConstrMeans[name_] = constrMean_;
    fConstrSigmas[name_] = constrSigma_;
    fConstrCorrParName[name_] = constrCorrParName_;
    fConstrCorrs[name_] = constrCorr_;

    AddParameter(name_, nom_, min_, max_, sigma_, nbins_, fdvalue_, label_);
  }
}
