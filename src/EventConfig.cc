#include <EventConfig.hh>
namespace antinufit
{
  double
  EventConfig::GetRate() const
  {
    return fRate;
  }

  void
  EventConfig::SetRate(double r_)
  {
    fRate = r_;
  }

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
  EventConfig::GetTexLabel() const
  {
    return fTexLabel;
  }

  void
  EventConfig::SetTexLabel(const std::string &s_)
  {
    fTexLabel = s_;
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

  std::string
  EventConfig::GetPdfPath() const
  {
    return fPdfPath;
  }

  void
  EventConfig::SetPdfPath(const std::string &s_)
  {
    fPdfPath = s_;
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

}
