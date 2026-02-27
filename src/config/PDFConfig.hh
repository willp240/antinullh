#pragma once

// OXO headers
#include <Exceptions.h>

// c++ headers
#include <vector>
#include <map>

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

    const std::vector<std::string> &GetDataBranchNames() const;
    void SetDataBranchNames(const std::vector<std::string> &);

    std::vector<std::string> GetBranchNames() const;
    std::vector<std::string> GetBranchNames(const int) const;

    const std::vector<std::string> &GetAxisNames() const;
    bool HasLLHBufferBins(const std::string &axisName) const;
    std::pair<int, int> GetLLHBufferBins(const std::string &axisName) const;
    void SetLLHBufferBins(const std::string &axisName, int low, int high);

  private:
    std::vector<std::string> fAxisNames;
    std::vector<std::string> fDataAxesNames;
    std::vector<std::string> fBranchNames;
    std::vector<std::string> fTexNames;
    std::vector<int> fBinCounts;
    std::vector<double> fMinima;
    std::vector<double> fMaxima;
    std::map<std::string, std::pair<int, int>> fLLHBufferBins = {{"energy", std::make_pair(8, 20)}};
  };
}
