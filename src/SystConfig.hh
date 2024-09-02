#ifndef __ANTINUFIT__SystConfig__
#define __ANTINUFIT__SystConfig__
#include <ParameterDict.h>
#include <string>
#include <vector>

namespace antinufit{
class SystConfig{
public:

  ParameterDict GetNominal() const;
  ParameterDict GetMinima() const;
  ParameterDict GetMaxima() const;
  ParameterDict GetMass() const;
  ParameterDict GetNBins() const;
  ParameterDict GetConstrMean() const;
  ParameterDict GetConstrSigma() const;
  ParameterDict GetSigma() const;
  std::map<std::string, std::string> GetGroup() const;
  std::map<std::string, std::string> GetType() const;
  std::map<std::string, std::string> GetObs() const;

  const std::string& GetName() const;
  void SetName(const std::string& name_);

  void AddParameter(const std::string& name_, double nom_, double min_, double max_, double mass_, double sigma_, int nbins_, const std::string& obs_, const std::string& type_, const std::string& group_);
  void AddParameter(const std::string& name_, double nom_, double min_, double max_, double mass, double sigma_, int nbins_, double constrMean_, double constrSigma_, const std::string& obs_, const std::string& type_, const std::string& group_);

private:
  std::string fName;
  ParameterDict fNominal;
  ParameterDict fMinima;  
  ParameterDict fMaxima;
  ParameterDict fMass;
  ParameterDict fConstrMean;
  ParameterDict fConstrSigma;
  ParameterDict fSigma;
  ParameterDict fNBins;
  std::map<std::string, std::string> fGroup;
  std::map<std::string, std::string> fType;
  std::map<std::string, std::string> fObs;
};
}
#endif
