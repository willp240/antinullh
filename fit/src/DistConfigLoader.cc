#include <DistConfigLoader.hh>
#include <DistConfig.hh>
#include <ConfigLoader.hh>
#include <map>
#include <algorithm>
#include <Exceptions.h>

namespace bbfit{

DistConfigLoader::DistConfigLoader(const std::string& filePath_){
    fPath = filePath_;    
}

DistConfig
DistConfigLoader::Load() const{
  ConfigLoader::Open(fPath);

  std::set<std::string> toLoad = ConfigLoader::ListSections();
  toLoad.erase("summary");

  std::vector<std::string> order;
  ConfigLoader::Load("summary", "build_order", order);

  std::string pdfDir;
  ConfigLoader::Load("summary", "pdf_dir", pdfDir);

  DistConfig retVal;

  double min;
  double max;
  std::string name;
  std::string branchName;
  std::string texName;
  int binCount;
  for(size_t i = 0; i < order.size(); i++){
    if(std::find(toLoad.begin(), toLoad.end(), order.at(i)) == toLoad.end())
      throw NotFoundError(Formatter() << "DistConfigLoader:: " << order.at(i)
			  << " is in the build order but has no section!");
    
    name = order.at(i);
    ConfigLoader::Load(name, "min", min);
    ConfigLoader::Load(name, "max", max);
    ConfigLoader::Load(name, "n_bins", binCount);
    ConfigLoader::Load(name, "branch_name", branchName);
    ConfigLoader::Load(name, "tex_name", texName);

    
    retVal.AddAxis(name, branchName, texName, binCount, min, max);
  }
  retVal.SetPDFDir(pdfDir);
  return retVal;
}

DistConfigLoader::~DistConfigLoader(){ 
  ConfigLoader::Close();
}

}
