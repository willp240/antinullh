#ifndef __BBFIT__DistConfig__
#define __BBFIT__DistConfig__
#include <string>
#include <vector>

namespace bbfit{
class DistConfig{
public:
  int GetAxisCount() const;
  void GetAxis(int index_, 
               std::string& name_, std::string& branchName_, std::string& texName_,
               int& binCount_, double& min_, double& max_) const;
  void AddAxis(const std::string& name_, const std::string& branchName_, 
               const std::string& texName_,
               int binCount_, double min_, double max_);

  const std::string& GetPDFDir() const;
  void SetPDFDir(const std::string&);

  const std::vector<std::string>& GetBranchNames() const;  

private:
  std::string fPDFDir;
  std::vector<std::string> fAxisNames;
  std::vector<std::string> fBranchNames;
  std::vector<std::string> fTexNames;
  std::vector<int> fBinCounts;
  std::vector<double> fMinima;
  std::vector<double> fMaxima;
};
}
#endif
