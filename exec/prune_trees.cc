/*
PruneTrees takes in ntuples and prunes out all the branches
we're not interested in so we only have nice lightweight files
to carry around in the fit.
It reads in input files specified in the event config and
literally just loops over events, filling new ntuples
with the quantities we want. These get written to wherever
was specified in the config file.
As the alpha-n's are simulated into the same files, we here
divide each alpha-n process into separate files. We use the
energy to determine which alpha-n process an event was
*/

// Antinu headers
#include <EventConfigLoader.hh>
#include <OscGridConfigLoader.hh>
#include <Utilities.hh>

// ROOT headers
#include <TChain.h>
#include <TFile.h>
#include <TNtuple.h>
#include <TVector3.h>

// c++ headers
#include <sys/stat.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <unordered_map>
#include <string>
#include <nlohmann/json.hpp>
#include <filesystem>

using namespace antinufit;

double precoil_cscatter_bound = 3.8;
double cscatter_oscatter_bound = 5.1;

void MakeDataSet(const std::vector<std::string> &filenames_,
                 const std::string &baseDir_,
                 const std::string &outFilename_,
                 std::unordered_map<std::string, int> &reactorNameIndex)
{

  // Output ntuple
  TFile outp(outFilename_.c_str(), "RECREATE");
  TNtuple nt("pruned", "", "energy:nu_energy:reactorIndex:alphaNClassifier");

  // Read the original data
  const std::string treeName = "output";
  for (size_t iFile = 0; iFile < filenames_.size(); iFile++)
  {
    std::string fileName = baseDir_ + "/" + filenames_.at(iFile);
    TChain chain(treeName.c_str());
    chain.Add(fileName.c_str());
    std::cout << fileName << "\t" << chain.GetEntries() << "  entries" << std::endl;

    int tenPercent = chain.GetEntries() / 10;

    Double_t energy;
    Double_t alphaNClassifier;
    TString *reactorName = NULL;
    Int_t reactorIndex;
    Double_t distance;
    Double_t neutrinoEnergy;

    chain.SetBranchAddress("energy", &energy);
    chain.SetBranchAddress("alphaNReactorIBD", &alphaNClassifier);
    chain.SetBranchAddress("parentKE1", &neutrinoEnergy);
    chain.SetBranchAddress("parentMeta1", &reactorName);

    // Read and write
    for (int i = 0; i < chain.GetEntries(); i++)
    {
      if (!(i % tenPercent))
      {
        std::cout << i << " / " << chain.GetEntries() << "\t ( " << 10 * i / tenPercent << " %)" << std::endl;
      }

      chain.GetEntry(i);

      std::string originReactorString(reactorName->Data());
      if (originReactorString != "")
      {
        reactorIndex = reactorNameIndex[originReactorString];
      }
      else
        reactorIndex = 999;

      // The alpha n particles are simulated at the same time into the same files. We'll split them now, by energy
      if (std::filesystem::path(outFilename_).filename().string() == "alphan_PRecoil.root" && energy > precoil_cscatter_bound)
        continue;
      if (std::filesystem::path(outFilename_).filename().string() == "alphan_CScatter.root" && (energy < precoil_cscatter_bound || energy > cscatter_oscatter_bound))
        continue;
      if (std::filesystem::path(outFilename_).filename().string() == "alphan_OExcited.root" && energy < cscatter_oscatter_bound)
        continue;

      outp.cd();
      nt.Fill(energy, neutrinoEnergy, reactorIndex, alphaNClassifier);
    }
  }
  outp.cd();
  nt.Write();
  return;
}

int main(int argc, char *argv[])
{
  if (argc != 3)
  {
    std::cout << "Usage: prune_trees <event_config_file> <osc_grid_config_file>" << std::endl;
    return 1;
  }

  // Load up the reactors json file
  std::string oscConfigFile(argv[2]);
  OscGridConfigLoader oscGridLoader(oscConfigFile);
  OscGridConfig oscGridConfig = oscGridLoader.Load();
  std::string reactorsJSONFile = oscGridConfig.GetReactorsJsonFile();
  std::unordered_map<std::string, int> reactorNameIndex = LoadNameIndexMap(reactorsJSONFile);

  // Create the results directory if it doesn't already exist
  std::string outDir;
  std::string eveConfigFile(argv[1]);
  ConfigLoader::Open(eveConfigFile);
  ConfigLoader::Load("summary", "pruned_ntup_dir", outDir);
  ConfigLoader::Close();
  EventConfigLoader eveLoader(eveConfigFile);
  typedef std::map<std::string, antinufit::EventConfig> EvMap;
  typedef std::vector<std::string> StringVec;
  typedef std::map<std::string, std::map<std::string, EventConfig>> DSMap;
  DSMap dsPDFMap = eveLoader.LoadActiveAndData();

  struct stat st = {0};
  if (stat(outDir.c_str(), &st) == -1)
  {
    mkdir(outDir.c_str(), 0700);
  }
  std::cout << "Made " << outDir << std::endl;

  for (DSMap::iterator dsIt = dsPDFMap.begin(); dsIt != dsPDFMap.end(); ++dsIt)
  {

    for (EvMap::iterator evIt = dsIt->second.begin(); evIt != dsIt->second.end(); ++evIt)
    {

      const std::string &name = evIt->first;
      std::cout << "Doing " << name << std::endl;
      const std::string &baseDir = evIt->second.GetNtupBaseDir();
      const StringVec &files = evIt->second.GetNtupFiles();
      const std::string &outName = evIt->second.GetPrunedPath();

      std::cout << "Writing from :" << std::endl;
      for (size_t i = 0; i < files.size(); i++)
        std::cout << "\t" << files.at(i) << std::endl;
      std::cout << "to " << outName << std::endl;

      MakeDataSet(files, baseDir, outName, reactorNameIndex);
    }
  }

  return 0;
}
