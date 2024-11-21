#ifndef __ANTINUFIT__SystConfig__
#define __ANTINUFIT__SystConfig__

// OXO headers
#include <ParameterDict.h>

// c++ headers
#include <vector>

namespace antinufit
{
  class SystConfig
  {
  public:
    std::map<std::string, std::string> GetGroup() const;
    std::map<std::string, std::string> GetType() const;
    std::map<std::string, std::vector<std::string>> GetDistObs() const;
    std::map<std::string, std::vector<std::string>> GetTransObs() const;
    std::map<std::string, std::string> GetParamNames() const;

    const std::string &GetName() const;
    void SetName(const std::string &name_);

    void AddParameter(const std::string &name_, const std::string &para_names_, const std::vector<std::string> &distObs_, const std::vector<std::string> &transObs_, const std::string &type_, const std::string &group_);

  private:
    std::string fName;
    std::map<std::string, std::string> fParamNames;
    std::map<std::string, std::string> fGroup;
    std::map<std::string, std::string> fType;
    std::map<std::string, std::vector<std::string> > fDistObs;
    std::map<std::string, std::vector<std::string> > fTransObs;
  };
}
#endif
