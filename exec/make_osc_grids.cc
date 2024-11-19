/*

*/

// Antinu headers
#include <OscGridConfigLoader.hh>
#include <OscGrid.hh>

// c++ headers
#include <sys/stat.h>
#include <string.h>
#include <map>
#include <tuple>
#include <string>
#include <fstream>
#include <nlohmann/json.hpp>

#include "TROOT.h"

#include <iostream>
#include <unistd.h>

using namespace antinufit;

int main(int argc, char *argv[])
{

  if (argc != 3)
  {
    std::cout << "Usage: make_osc_grids oscgrid_config index_to_write" << std::endl;
    return 1;
  }

  std::string oscgridConfigFile(argv[1]);
  int index = std::stoi(argv[2]);
  OscGridConfigLoader oscGridLoader(oscgridConfigFile);
  OscGridConfig oscGridConfig = oscGridLoader.Load();

  std::string outfilename = oscGridConfig.GetFilename();
  std::string reactorsjsonfile = oscGridConfig.GetReactorsJsonFile();

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
  std::unordered_map<int, double> indexDistance = LoadIndexDistanceMap(reactorsjsonfile);

  std::string oscGridFileName = outfilename + "_" + std::to_string(index) + ".root";
  std::cout << "Making grid " << oscGridFileName << std::endl;
  double distance = indexDistance[index];
  std::unique_ptr<OscGrid> oscGrid = std::make_unique<OscGrid>(oscGridFileName, distance, minE, maxE, numValsE, minDm21sq, maxDm21sq, numValsDm21sq, minSsqth12, maxSsqth12, numValsSsqth12);
  oscGrid->CalcGrid();
  oscGrid->Write();

return 0;
}