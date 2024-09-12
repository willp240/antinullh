#ifndef __ANTINUFIT__FitConfig__
#define __ANTINUFIT__FitConfig__
#include <ParameterDict.h>
#include <string>
#include <set>

namespace antinufit
{
  class FitConfig
  {
  public:
    ParameterDict GetMinima() const;
    ParameterDict GetMaxima() const;
    ParameterDict GetSigmas() const;
    ParameterDict GetNBins() const;
    ParameterDict GetNominals() const;

    ParameterDict GetConstrMeans() const;
    ParameterDict GetConstrSigmas() const;

    int GetIterations() const;
    void SetIterations(int);
    int GetHMCIterations() const;
    void SetHMCIterations(int);

    int GetBurnIn() const;
    void SetBurnIn(int);
    int GetHMCBurnIn() const;
    void SetHMCBurnIn(int);

    void AddParameter(const std::string &name_, double mean_, double min_, double max_, double sigma_, int nbins_);
    void AddParameter(const std::string &name_, double mean_, double min_, double max_, double sigma_, int nbins_,
                      double constrMean_, double constrSigma_);

    std::set<std::string> GetParamNames() const;

    const std::string &GetOutDir() const;
    void SetOutDir(const std::string &);

    int GetNSteps() const;
    void SetNSteps(int);

    double GetEpsilon() const;
    void SetEpsilon(double);

    double GetSigmaScale() const;
    void SetSigmaScale(double);

    bool GetBeestonBarlow() const;
    void SetBeestonBarlow(bool);

    bool GetAsimov() const;
    void SetAsimov(bool);

    std::string GetDatafile() const;
    void SetDatafile(std::string);

    double GetLivetime() const;
    void SetLivetime(double);

  private:
    std::string fOutDir;
    ParameterDict fConstrMeans;
    ParameterDict fConstrSigmas;
    ParameterDict fNominals;
    ParameterDict fMinima;
    ParameterDict fMaxima;
    ParameterDict fSigmas;
    ParameterDict fNbins;
    int fIterations;
    int fBurnIn;
    int fHMCIterations;
    int fHMCBurnIn;
    int fNsteps;
    bool fBeestonBarlow;
    double fEpsilon;
    double fSigmaScale;
    bool fAsimov;
    double fLivetime;
    std::string fDatafile;
  };
}
#endif
