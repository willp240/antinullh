#include <CutConfigLoader.hh>
#include <CutConfig.hh>
#include <ConfigLoader.hh>
#include <map>
#include <algorithm>

namespace bbfit{

CutConfigLoader::CutConfigLoader(const std::string& filePath_){
    fPath = filePath_;    
}

CutConfig
CutConfigLoader::LoadOne(const std::string& name_) const{
  ConfigLoader::Open(fPath);

  double value;
  double value2;
  std::string type;
  std::string obs;

  ConfigLoader::Load(name_, "value", value);  
  ConfigLoader::Load(name_, "type", type);
  ConfigLoader::Load(name_, "obs", obs);
  
  try{
      ConfigLoader::Load(name_, "value2", value2);
  }
  catch(const ConfigFieldMissing& ){
      if(type != "bool")
          throw;
      
  }
  
  CutConfig retVal;
  retVal.SetName(name_);
  retVal.SetValue(value);
  retVal.SetValue2(value2);
  retVal.SetType(type);
  retVal.SetObs(obs);
  return retVal;
}

std::vector<CutConfig>
CutConfigLoader::LoadActive() const{
  typedef std::vector<std::string> StringVec;
  typedef std::vector<CutConfig> CutConfigVec;
  ConfigLoader::Open(fPath);

  StringVec toLoad;
  ConfigLoader::Load("summary", "order", toLoad);

  CutConfigVec evVec;
  for(size_t i = 0; i < toLoad.size(); i++)
    evVec.push_back(LoadOne(toLoad.at(i)));

  return evVec;
}

CutConfigLoader::~CutConfigLoader(){ 
  ConfigLoader::Close();
}

}
