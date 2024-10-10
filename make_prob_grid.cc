/*

*/

// Antinu headers
#include <OscGridConfigLoader.hh>

// c++ headers
#include <sys/stat.h>
#include <string.h>
#include <iostream>
#include <map>
#include <tuple>
#include <string>
#include <fstream>
#include <nlohmann/json.hpp>

typedef std::map<std::string, std::string> StringMap;

int main(int argc, char *argv[])
{
  if (argc != 2)
  {
    std::cout << "Usage: make_prob_grid oscgrid_config" << std::endl;
    return 1;
  }

  std::string oscgridConfigFile(argv[1]);
  OscGridConfigLoader oscgridLoader(oscgridConfigFile);
  OscGridConfig oscGridConfig = oscGridLoader.Load();

  std::string outfilename = oscGridConfig.GetFilename();
  double distance = oscGridConfig.GetDistance();
  double minE = oscGridConfig.GetMinE();
  double maxE = oscGridConfig.GetMaxE();
  int numValsE = oscGridConfig.GetNumValsE();
  double minDm21sq = oscGridConfig.GetMinDm21sq();
  double maxDm21sq = oscGridConfig.GetMaxDm21sq();
  int numValsDm21sq = oscGridConfig.GetNumValsDm21sq();
  double minSsqth12 = oscGridConfig.GetMinSsqth12();
  double maxSsqth12 = oscGridConfig.GetMaxSsqth12();
  int numValsSsqth12 = oscGridConfig.GetNumValsSsqth12();

  // First read the reactor distance info
  std::unordered_map<int, double> indexDistance = LoadIndexDistanceMap("reactors.json");

  for (std::map<int, double>::iterator it = indexDistance.begin(); it != indexDistance.end(); ++it)
  {
    std::string oscGridFileName = outfilename + "_" + std::to_string(it->first) + ".csv";
    OscGrid oscGrid(outfilename, distance, minE, maxE, numValsE, minDm21sq, numValsDm21sq, minSsqth12, maxSsqth12, numValsSsqth12);
    oscGrid.CalcGrid();
    oscGrid.Write();
  }

}