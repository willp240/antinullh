/*
Takes in an OscGridConfig (for the JSON filename) and loops over reactor cores
in the REACTOR RATDB table and calculates the distance. Then assigns a unique
identifying int for the core, and writes each combination of reactor name,
distance, and int to the JSON file specified in the Osc Grid config.
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

// AntinuLLH headers
#include <OscGridConfigLoader.hh>
#include <Utilities.hh>

using namespace antinufit;

int main(int argc, char *argv[])
{
  if (argc != 2)
  {
    std::cout << "Usage: make_reactor_json osc_grid_config_file" << std::endl;
    return 1;
  }

  // Load up the reactors json file
  std::string oscConfigFile(argv[2]);
  OscGridConfigLoader oscGridLoader(oscConfigFile);
  OscGridConfig oscGridConfig = oscGridLoader.Load();
  std::string reactorsJSONFile = oscGridConfig.GetReactorsJsonFile();

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
    std::vector<Double_t> latitude = linkdb->GetDArray("latitude");
    std::vector<Double_t> longitute = linkdb->GetDArray("longitude");
    std::vector<Double_t> altitude = linkdb->GetDArray("altitude");
    for (int iCore = 0; iCore < numCores; iCore++)
    {
      std::string reactCoreName = it->first + " " + std::to_string(iCore);
      std::cout << reactCoreName << std::endl;

      const double baseline = GetReactorDistanceLLA(longitute[iCore], latitude[iCore], altitude[iCore]);
      reactorIndex[index] = std::make_tuple(reactCoreName, baseline);
      index++;
    }
  }

  nlohmann::json j;

  for (const auto &[key, value] : reactorIndex)
  {
    j[std::to_string(key)] = {std::get<0>(value), std::get<1>(value)};
  }

  // Save to file
  std::ofstream file(reactorsJSONFile);
  file << j.dump(4); // Pretty-print with 4 spaces
  file.close();

  return 0;
}
