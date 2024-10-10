/*

*/

// c++ headers
#include <sys/stat.h>
#include <string.h>
#include <iostream>
#include <map>
#include <tuple>
#include <string>
#include <fstream>
#include <nlohmann/json.hpp>

// RAT headers
#include <RAT/DB.hh>

typedef std::map<std::string, std::string> StringMap;

TVector3 LLAtoECEF(double longitude, double latitude, double altitude)
{
  // reference http://www.mathworks.co.uk/help/aeroblks/llatoecefposition.html
  static double toRad = TMath::Pi() / 180.;
  static double Earthradius = 6378137.0; // Radius of the Earth (in meters)
  static double f = 1. / 298.257223563;  // Flattening factor WGS84 Model
  static double L, rs, x, y, z;
  L = atan(pow((1. - f), 2) * tan(latitude * toRad)) * 180. / TMath::Pi();
  rs = sqrt(pow(Earthradius, 2) / (1. + (1. / pow((1. - f), 2) - 1.) * pow(sin(L * toRad), 2)));
  x = (rs * cos(L * toRad) * cos(longitude * toRad) + altitude * cos(latitude * toRad) * cos(longitude * toRad)) / 1000; // in km
  y = (rs * cos(L * toRad) * sin(longitude * toRad) + altitude * cos(latitude * toRad) * sin(longitude * toRad)) / 1000; // in km
  z = (rs * sin(L * toRad) + altitude * sin(latitude * toRad)) / 1000;                                                   // in km

  TVector3 ECEF = TVector3(x, y, z);

  return ECEF;
}


double GetReactorDistanceLLA(const double &longitude, const double &latitude, const double &altitude)
{
  const TVector3 SNO_ECEF_coord_ = TVector3(672.87, -4347.18, 4600.51);
  double dist = (LLAtoECEF(longitude, latitude, altitude) - SNO_ECEF_coord_).Mag();
  return dist;
}


std::vector<std::string> SplitString(std::string str, bool Osc){
    std::istringstream buf(str);
    std::istream_iterator<std::string> beg(buf), end;
    std::vector<std::string> tokens(beg, end); //each word of string now in vector
    std::vector<std::string> info;
    std::string dummy = "";

    for(int i=0;i<tokens.size();i++){ //combine back to reactor name, core number
        if(i==0){
            dummy = tokens.at(i);
            if(i!=tokens.size()-2){
                dummy += " ";
            }
        }
        else if(i!=tokens.size()-1){
            dummy += tokens.at(i);
            if(i!=tokens.size()-2){
                dummy += " ";
            }
        }
        else{
            info.push_back(dummy);
            info.push_back(tokens.at(i));
        }
    }

    return info;
}


int main(int argc, char *argv[])
{
  if (argc != 1)
  {
    std::cout << "Usage: make_reactor_json" << std::endl;
    return 1;
  }

  std::map<int, std::tuple<std::string, double>> reactorIndex;

  RAT::DBLinkPtr linkdb;
  RAT::DB *db = RAT::DB::Get();
  db->LoadDefaults();

  RAT::DBLinkGroup grp = db->GetLinkGroup("REACTOR");
  RAT::DBLinkGroup::iterator it;
  int index = 0;
  for (it = grp.begin(); it != grp.end(); ++it)
  {
    linkdb = RAT::DB::Get()->GetLink("REACTOR", it->first);
    Double_t numCores = linkdb->GetD("no_cores");
    std::vector<Double_t> latitude  = linkdb->GetDArray("latitude");
    std::vector<Double_t> longitute = linkdb->GetDArray("longitude");
    std::vector<Double_t> altitude = linkdb->GetDArray("altitude");
    for (int iCore = 0; iCore < numCores; iCore++){
      std::string reactCoreName = it->first + " " + std::to_string(iCore);
      std::cout << reactCoreName << std::endl;

      const double baseline = GetReactorDistanceLLA(longitute[iCore], latitude[iCore], altitude[iCore]);
      reactorIndex[index] = std::make_tuple(reactCoreName, baseline);
      index++;

    }
  }

  nlohmann::json j;

  for (const auto& [key, value] : reactorIndex) {
      j[std::to_string(key)] = { std::get<0>(value), std::get<1>(value) };
  }

  // Save to file
  std::ofstream file("reactors.json");
  file << j.dump(4);  // Pretty-print with 4 spaces
  file.close();

return 0;
}
