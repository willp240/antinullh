#include <SystConfigLoader.hh>
#include <SystConfig.hh>
#include <ConfigLoader.hh>
#include <map>
#include <algorithm>

namespace antinufit{

SystConfigLoader::SystConfigLoader(const std::string& filePath_){
    fPath = filePath_;    
}

SystConfig
SystConfigLoader::LoadActive() const{
  SystConfig ret;
  ConfigLoader::Open(fPath);
 
  typedef std::set<std::string> StringSet;
  StringSet toLoad;
  ConfigLoader::Load("summary", "active", toLoad);

  std::string name;
  double      nom;
  double      min;
  double      max;
  double      mass;
  double      sigma;
  double      constrMean;
  double      constrSigma;
  int         nbins;
  double      std_nom;
  double      std_min;
  double      std_max;
  double      std_mass;
  double      std_sigma;
  int         std_nbins;
  std::string obs;
  std::string type;

  if(std::find(toLoad.begin(), toLoad.end(), "all") != toLoad.end()){
    toLoad = ConfigLoader::ListSections();
    toLoad.erase("summary");
  }

  for(StringSet::iterator it = toLoad.begin(); it != toLoad.end(); ++it){
    name = *it;
    ConfigLoader::Load(name, "nominal", nom);
    ConfigLoader::Load(name, "minima", min);
    ConfigLoader::Load(name, "maxima", max);
    ConfigLoader::Load(name, "mass", mass);
    ConfigLoader::Load(name, "sigma", sigma);
    ConfigLoader::Load(name, "nbins", nbins);
    ConfigLoader::Load(name, "obs", obs);
    ConfigLoader::Load(name, "type", type);

    if(type!="convolution"){ 
      try{
	ConfigLoader::Load(name, "constraint_mean", constrMean);
	ConfigLoader::Load(name, "constraint_sigma", constrSigma);
	ret.AddParameter(name, nom, min, max, mass, sigma, nbins, constrMean, constrSigma, obs, type);
      }
      catch(const ConfigFieldMissing& e_){
	ret.AddParameter(name, nom, min, max, mass, sigma, nbins, obs, type);
      }
    }
    else{
      ConfigLoader::Load(name, "nominal_stddev", std_nom);
      ConfigLoader::Load(name, "minima_stddev", std_min);
      ConfigLoader::Load(name, "maxima_stddev", std_max);
      ConfigLoader::Load(name, "mass_stddev", std_mass);
      ConfigLoader::Load(name, "nbins_stddev", std_nbins);
      try{
	ConfigLoader::Load(name, "constraint_mean", constrMean);
	ConfigLoader::Load(name, "constraint_sigma", constrSigma);
	ret.AddParameter(name, nom, min, max, mass, sigma, nbins, constrMean, constrSigma, obs, type, std_nom, std_min, std_max, std_mass, std_sigma, std_nbins);
      }
      catch(const ConfigFieldMissing& e_){
        ret.AddParameter(name, nom, min, max, mass, sigma, nbins, obs, type, std_nom, std_min, std_max, std_mass, std_sigma, std_nbins);
      }
    }
  }
  return ret;
}

SystConfigLoader::~SystConfigLoader(){ 
  ConfigLoader::Close();
}

}
