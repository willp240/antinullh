#include <SystConfig.hh>

namespace antinufit
{

  std::map<std::string, std::vector<std::string> >
  SystConfig::GetDistObs() const
  {
    return fDistObs;
  }

  std::map<std::string, std::vector<std::string> >
  SystConfig::GetTransObs() const
  {
    return fTransObs;
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

  const std::string &
  SystConfig::GetName() const
  {
    return fName;
  }

  void
  SystConfig::AddParameter(const std::string &name_, const std::string &param_names_, const std::vector<std::string> &dist_obs_, const std::vector<std::string> &trans_obs_, const std::string &type_, const std::string &group_)
  {
    fDistObs[name_] = dist_obs_;
    fTransObs[name_] = trans_obs_;
    fType[name_] = type_;
    fGroup[name_] = group_;
    fParamNames[name_] = param_names_;
  }
}
