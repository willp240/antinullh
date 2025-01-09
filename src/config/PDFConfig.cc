#include <PDFConfig.hh>

namespace antinufit
{
  int PDFConfig::GetAxisCount() const
  {
    return fAxisNames.size();
  }

  int PDFConfig::GetDataAxisCount() const
  {
    return fDataAxesNames.size();
  }

  void
  PDFConfig::GetAxis(int index_, std::string &name_,
                      std::string &branchName_,
                      std::string &texName_,
                      int &binCount_, double &min_, double &max_) const
  {
    try
    {
      name_ = fAxisNames.at(index_);
      branchName_ = fBranchNames.at(index_);
      texName_ = fTexNames.at(index_);
      binCount_ = fBinCounts.at(index_);
      min_ = fMinima.at(index_);
      max_ = fMaxima.at(index_);
    }
    catch (const std::out_of_range &e_)
    {
      throw NotFoundError(Formatter() << "PDFConfig::No data for axis " << index_);
    }
  }

  void
  PDFConfig::AddAxis(const std::string &name_, const std::string &branchName_,
                      const std::string &texName_,
                      int binCount_, double min_, double max_)
  {
    fAxisNames.push_back(name_);
    fTexNames.push_back(texName_);
    fBranchNames.push_back(branchName_);
    fBinCounts.push_back(binCount_);
    fMinima.push_back(min_);
    fMaxima.push_back(max_);
  }

  const std::string &
  PDFConfig::GetPDFDir() const
  {
    return fPDFDir;
  }

  void
  PDFConfig::SetPDFDir(const std::string &s_)
  {
    fPDFDir = s_;
  }

  std::vector<std::string>
  PDFConfig::GetBranchNames() const
  {
    int numBranches = fBranchNames.size();
    return GetBranchNames(numBranches);
  }

  std::vector<std::string>
  PDFConfig::GetBranchNames( const int numBranches_) const
  {
    std::vector<std::string> bnames(fBranchNames.begin(), fBranchNames.begin() + numBranches_);
    return bnames;
  }

  void
  PDFConfig::SetDataBranchNames(const std::vector<std::string> &s_)
  {
    fDataAxesNames = s_;
  }

  const std::vector<std::string> &
  PDFConfig::GetDataBranchNames() const
  {
    return fDataAxesNames;
  }

}
