/*

*/

// Antinu headers
#include <OscGridConfigLoader.hh>
#include <OscGrid.hh>

// c++ headers
#include <sys/stat.h>
#include <string.h>
#include <iostream>
#include <map>
#include <tuple>
#include <string>
#include <fstream>
#include <nlohmann/json.hpp>
#include <unistd.h>

#include <TFile.h>
#include <TH3D.h>

// typedef std::map<std::string, std::string> StringMap;

using namespace antinufit;

int main(int argc, char *argv[])
{
  std::cout << "pre args" << std::endl;

  if (argc != 2)
  {
    std::cout << "Usage: make_osc_grids oscgrid_config" << std::endl;
    return 1;
  }

  std::string oscgridConfigFile(argv[1]);
  OscGridConfigLoader oscGridLoader(oscgridConfigFile);
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
  std::unordered_map<int, double> indexDistance = LoadIndexDistanceMap("reactors_bruce1.json");

  for (std::unordered_map<int, double>::iterator it = indexDistance.begin(); it != indexDistance.end(); ++it)
  {
    std::string oscGridFileName = outfilename + "_" + std::to_string(it->first) + ".csv";
    std::cout << "Making grid " << oscGridFileName << std::endl;
    OscGrid *oscGrid = new OscGrid(oscGridFileName, it->second, minE, maxE, numValsE, minDm21sq, maxDm21sq, numValsDm21sq, minSsqth12, maxSsqth12, numValsSsqth12);
    oscGrid->CalcGrid();
    oscGrid->Write();
    std::cout << oscGrid->Evaluate(5, 7.5E-5, 0.5) << std::endl;

    TH3D *h3 = oscGrid->MakeHist();
    TFile *file = new TFile("hist140x50x50.root", "RECREATE");

    h3->Write();
    file->Close();
    delete oscGrid;
  }
}