#ifndef __BBFIT__EventConfig__
#define __BBFIT__EventConfig__
#include <string>
#include <vector>

namespace bbfit{
class EventConfig{
public:
  EventConfig(): fRate(-1), fNgenerated(0) {}

  double GetRate() const;
  void   SetRate(double);

  unsigned long GetNGenerated() const;
  void SetNGenerated(unsigned long);

  std::string GetName() const;
  void SetName(const std::string&);

  const std::vector<std::string>& GetNtupFiles() const;
  void SetNtupFiles(const std::vector<std::string>&);

  std::string GetTexLabel() const;
  void SetTexLabel(const std::string&);

  std::string GetLoadingScaling() const;
  void SetLoadingScaling(const std::string&);

  std::string GetNtupBaseDir() const;
  void SetNtupBaseDir(const std::string&);

  std::string GetPrunedPath() const;
  void SetPrunedPath(const std::string&);

  std::string GetSplitFakePath() const;
  void SetSplitFakePath(const std::string&);

	std::string GetSplitPdfPath() const;
  void SetSplitPdfPath(const std::string&);

  bool GetRandomSplit() const;
  void SetRandomSplit(bool);

private:
  double fRate;
  unsigned long  fNgenerated;
  std::vector<std::string> fNtupFiles;
  std::string fNtupBaseDir; // the originals
  std::string fPrunedPath;   // the pruned ouput
  std::string fSplitFakePath;   // after the fake data split
	std::string fSplitPdfPath;   // after the fake data split
  std::string fTexLabel;
  std::string fName;
  std::string fLoadingScaling;
  bool fRandomSplit;
};
}

#endif
