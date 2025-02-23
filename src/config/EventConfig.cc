#include <EventConfig.hh>

namespace antinufit
{

  const std::vector<std::string> &
  EventConfig::GetNtupFiles() const
  {
    return fNtupFiles;
  }

  void
  EventConfig::SetNtupFiles(const std::vector<std::string> &s_)
  {
    fNtupFiles = s_;
  }

  std::string
  EventConfig::GetName() const
  {
    return fName;
  }

  void
  EventConfig::SetName(const std::string &s_)
  {
    fName = s_;
  }

  std::string
  EventConfig::GetNtupBaseDir() const
  {
    return fNtupBaseDir;
  }

  void
  EventConfig::SetNtupBaseDir(const std::string &s_)
  {
    fNtupBaseDir = s_;
  }

  std::string
  EventConfig::GetPrunedPath() const
  {
    return fPrunedPath;
  }

  void
  EventConfig::SetPrunedPath(const std::string &s_)
  {
    fPrunedPath = s_;
  }

  std::vector<std::string>
  EventConfig::GetGroup() const
  {
    return fGroup;
  }

  void
  EventConfig::SetGroup(const std::vector<std::string> &s_)
  {
    fGroup = s_;
  }

  int
  EventConfig::GetNumDimensions() const
  {
    return fNumDimensions;
  }

  void
  EventConfig::SetNumDimensions(const int &i_)
  {
    fNumDimensions = i_;
  }

}
