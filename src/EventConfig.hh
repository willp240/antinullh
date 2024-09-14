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
    EventConfig() : fRate(-1) {}

    double GetRate() const;
    void SetRate(double);

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

    std::string GetPdfPath() const;
    void SetPdfPath(const std::string &);

    std::vector<std::string> GetGroup() const;
    void SetGroup(const std::vector<std::string> &);

  private:
    double fRate;
    std::vector<std::string> fNtupFiles;
    std::string fNtupBaseDir; // the originals
    std::string fPrunedPath;  // the pruned ouput
    std::string fPdfPath;
    std::string fTexLabel;
    std::string fName;
    std::vector<std::string> fGroup;
  };
}

#endif
