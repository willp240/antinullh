#ifndef __ANTINUFIT__PDFConfig__
#define __ANTINUFIT__PDFConfig__

// OXO headers
#include <Exceptions.h>

// c++ headers
#include <vector>

namespace antinufit
{
  class PDFConfig
  {
  public:
    int GetAxisCount() const;
    int GetDataAxisCount() const;
    void GetAxis(int index_,
                 std::string &name_, std::string &branchName_, std::string &texName_,
                 int &binCount_, double &min_, double &max_) const;
    void AddAxis(const std::string &name_, const std::string &branchName_,
                 const std::string &texName_,
                 int binCount_, double min_, double max_);

    const std::string &GetPDFDir() const;
    void SetPDFDir(const std::string &);

    const std::vector<std::string> &GetDataBranchNames() const;
    void SetDataBranchNames(const std::vector<std::string> &);

    std::vector<std::string> GetBranchNames( ) const;
    std::vector<std::string> GetBranchNames( const int ) const;

  private:
    std::string fPDFDir;
    std::vector<std::string> fAxisNames;
    std::vector<std::string> fDataAxesNames;
    std::vector<std::string> fBranchNames;
    std::vector<std::string> fTexNames;
    std::vector<int> fBinCounts;
    std::vector<double> fMinima;
    std::vector<double> fMaxima;
  };
}
#endif
