#ifndef __ANTINUFIT__FitConfig__
#define __ANTINUFIT__FitConfig__
#include <ParameterDict.h>
#include <string>
#include <set>

namespace antinufit{
  class FitConfig{
  public:
    ParameterDict GetMinima() const;
    ParameterDict GetMaxima() const;
    ParameterDict GetSigmas() const;
    ParameterDict GetNBins() const;
  
    ParameterDict GetConstrMeans() const;
    ParameterDict GetConstrSigmas() const;

    int  GetIterations() const;
    void SetIterations(int);
    int  GetHMCIterations() const;
    void SetHMCIterations(int);

    int  GetBurnIn() const;
    void SetBurnIn(int);
    int  GetHMCBurnIn() const;
    void SetHMCBurnIn(int);

    void AddParameter(const std::string& name_, double min_, double max_, double sigma_, int nbins_);
    void AddParameter(const std::string& name_, double min_, double max_, double sigma_, int nbins_, 
		      double constrMean_, double constrSigma_);

    std::set<std::string> GetParamNames() const;
  
    const std::string& GetOutDir() const;
    void  SetOutDir(const std::string&);

    int  GetNSteps() const;
    void SetNSteps(int);

    double GetEpsilon() const;
    void   SetEpsilon(double);

    double GetSigmaScale() const;
    void   SetSigmaScale(double);

    bool   GetBeestonBarlow() const;
    void   SetBeestonBarlow(bool);

  private:
    std::string   fOutDir;
    ParameterDict fConstrMeans;
    ParameterDict fConstrSigmas;;
    ParameterDict fMinima;
    ParameterDict fMaxima;
    ParameterDict fSigmas;
    ParameterDict fNbins;
    int       fIterations;
    int       fBurnIn;
    int       fHMCIterations;
    int       fHMCBurnIn;
    int       fNsteps;
    bool      fBeestonBarlow;
    double    fEpsilon;
    double    fSigmaScale;
    
  };
}
#endif

