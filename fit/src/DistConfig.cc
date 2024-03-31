#include <DistConfig.hh>
#include <Exceptions.h>


namespace bbfit{
int
DistConfig::GetAxisCount() const{
  return fAxisNames.size();
}

void
DistConfig::GetAxis(int index_, std::string& name_, 
		    std::string& branchName_, 
		    std::string& texName_, 
		    int& binCount_, double& min_, double& max_) const{
  try{
    name_ = fAxisNames.at(index_);
    branchName_ = fBranchNames.at(index_);
    texName_ = fTexNames.at(index_);
    binCount_ = fBinCounts.at(index_);
    min_ = fMinima.at(index_);
    max_ = fMaxima.at(index_);
  }
  catch(const std::out_of_range& e_){
    throw NotFoundError(Formatter() << "DistConfig::No data for axis " << index_);
  }
}

void
DistConfig::AddAxis(const std::string& name_, const std::string& branchName_, 
		   const std::string& texName_, 
		   int binCount_, double min_, double max_){
  fAxisNames.push_back(name_);
  fTexNames.push_back(texName_);
  fBranchNames.push_back(branchName_);
  fBinCounts.push_back(binCount_);
  fMinima.push_back(min_);
  fMaxima.push_back(max_);  
}

const std::string&
DistConfig::GetPDFDir() const{
  return fPDFDir;
}

void
DistConfig::SetPDFDir(const std::string& s_){
  fPDFDir = s_;
}

const std::vector<std::string>&
DistConfig::GetBranchNames() const {
  return fBranchNames;
}

}
