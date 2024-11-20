#ifndef __ANTINUFIT__EventConfig__
#define __ANTINUFIT__EventConfig__

// c++ headers
#include <string>
#include <vector>

namespace antinufit
{
  class EventConfig
  {
  public:
    EventConfig() {}

    std::string GetName() const;
    void SetName(const std::string &);

    const std::vector<std::string> &GetNtupFiles() const;
    void SetNtupFiles(const std::vector<std::string> &);

    std::string GetTexLabel() const;
    void SetTexLabel(const std::string &);

    std::string GetNtupBaseDir() const;
    void SetNtupBaseDir(const std::string &);

    std::string GetPrunedPath() const;
    void SetPrunedPath(const std::string &);

    std::vector<std::string> GetGroup() const;
    void SetGroup(const std::vector<std::string> &);

    int GetNumDimensions() const;
    void SetNumDimensions(const int &);

  private:
    std::vector<std::string> fNtupFiles;
    std::string fNtupBaseDir; // The originals
    std::string fPrunedPath;  // The pruned ouput
    std::string fTexLabel;
    std::string fName;
    std::vector<std::string> fGroup;
    int fNumDimensions;
  };
}

#endif
