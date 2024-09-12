#include <SystConfig.hh>
#include <ContainerTools.hpp>

namespace antinufit
{

  std::map<std::string, std::string>
  SystConfig::GetObs() const
  {
    return fObs;
  }

  std::map<std::string, std::string>
  SystConfig::GetType() const
  {
    return fType;
  }

  std::map<std::string, std::string>
  SystConfig::GetGroup() const
  {
    return fGroup;
  }

  std::map<std::string, std::string>
  SystConfig::GetParamNames() const
  {
    return fParamNames;
  }

  std::map<std::string, std::string>
  SystConfig::GetFunctionNames() const
  {
    return fFunctionNames;
  }

  const std::string &
  SystConfig::GetName() const
  {
    return fName;
  }

  void
  SystConfig::AddParameter(const std::string &name_, const std::string &param_names_, const std::string &obs_, const std::string &type_, const std::string &group_, const std::string &function_)
  {
    fObs[name_] = obs_;
    fType[name_] = type_;
    fGroup[name_] = group_;
    fParamNames[name_] = param_names_;
    fFunctionNames[name_] = function_;
  }
}
