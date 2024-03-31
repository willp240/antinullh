#include <EventConfig.hh>
namespace bbfit{
double
EventConfig::GetRate() const{
  return fRate;
}

void
EventConfig::SetRate(double r_){
  fRate = r_;
}

unsigned long
EventConfig::GetNGenerated() const{
  return fNgenerated;
}

void
EventConfig::SetNGenerated(unsigned long n_){
  fNgenerated = n_;

}

const std::vector<std::string>&
EventConfig::GetNtupFiles() const{
  return fNtupFiles;
}

void
EventConfig::SetNtupFiles(const std::vector<std::string>& s_){
  fNtupFiles = s_;
}

std::string
EventConfig::GetTexLabel() const{
  return fTexLabel;
}

void
EventConfig::SetTexLabel(const std::string& s_){
  fTexLabel = s_;
}

std::string
EventConfig::GetLoadingScaling() const{
  return fLoadingScaling;
}

void
EventConfig::SetLoadingScaling(const std::string& s_){
  fLoadingScaling = s_;
}

std::string
EventConfig::GetName() const{
  return fName;
}

void
EventConfig::SetName(const std::string& s_){
  fName = s_;
}

std::string
EventConfig::GetNtupBaseDir() const{
  return fNtupBaseDir;
}

void
EventConfig::SetNtupBaseDir(const std::string& s_){
  fNtupBaseDir = s_;
}

std::string
EventConfig::GetPrunedPath() const{
  return fPrunedPath;
}

void
EventConfig::SetPrunedPath(const std::string& s_){
  fPrunedPath = s_;
}

std::string
EventConfig::GetSplitFakePath() const{
  return fSplitFakePath;
}

void
EventConfig::SetSplitFakePath(const std::string& s_){
  fSplitFakePath = s_;
}
	
std::string
EventConfig::GetSplitPdfPath() const{
  return fSplitPdfPath;
}

void
EventConfig::SetSplitPdfPath(const std::string& s_){
  fSplitPdfPath = s_;
}

bool
EventConfig::GetRandomSplit() const{
  return fRandomSplit;
}

void
EventConfig::SetRandomSplit(bool b_){
  fRandomSplit = b_;
}
}
