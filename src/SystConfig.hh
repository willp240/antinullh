#ifndef __ANTINUFIT__SystConfig__
#define __ANTINUFIT__SystConfig__

// OXO headers
#include <ParameterDict.h>

namespace antinufit
{
  class SystConfig
  {
  public:
    std::map<std::string, std::string> GetGroup() const;
    std::map<std::string, std::string> GetType() const;
    std::map<std::string, std::string> GetObs() const;
    std::map<std::string, std::string> GetParamNames() const;
    std::map<std::string, std::string> GetFunctionNames() const;

    const std::string &GetName() const;
    void SetName(const std::string &name_);

    void AddParameter(const std::string &name_, const std::string &para_names_, const std::string &obs_, const std::string &type_, const std::string &group_, const std::string &function_);

  private:
    std::string fName;
    std::map<std::string, std::string> fParamNames;
    std::map<std::string, std::string> fGroup;
    std::map<std::string, std::string> fType;
    std::map<std::string, std::string> fObs;
    std::map<std::string, std::string> fFunctionNames;
  };
}
#endif
