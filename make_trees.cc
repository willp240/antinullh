/*
MakeTrees takes in ntuples and prunes out all the branches
we're not interested in so we only have nice lightweight files
to carry around in the fit.
It reads in input files specified in the event config and
literally just loops over events, filling new ntuples
with the quantities we want. These get written to wherever
was specified in the config file.
*/

#include <TChain.h>
#include <TFile.h>
#include <TNtuple.h>
#include <string>
#include <vector>
#include <iostream>
#include <cmath>
#include <sstream>
#include <map>
#include <EventConfigLoader.hh>
#include <ConfigLoader.hh>
#include <sys/stat.h>

typedef std::map<std::string, std::string> StringMap;

void MakeDataSet(const std::vector<std::string> &filenames_,
                 const std::string &baseDir_,
                 const std::string &treeName_,
                 const std::string &outFilename_)
{

  // output ntuple
  TFile outp(outFilename_.c_str(), "RECREATE");
  TNtuple nt("pruned", "", "energy:alphaNReactorIBD");

  // read the original data
  for (size_t iFile = 0; iFile < filenames_.size(); iFile++)
  {
    std::string fileName = baseDir_ + "/" + filenames_.at(iFile);
    TChain chain(treeName_.c_str());
    chain.Add(fileName.c_str());
    std::cout << fileName << "\t" << chain.GetEntries() << "  entries" << std::endl;

    int tenPercent = chain.GetEntries() / 10;

    Double_t e;
    Double_t c;

    chain.SetBranchAddress("energy", &e);
    chain.SetBranchAddress("alphaNReactorIBD", &c);

    // read and write

    for (int i = 0; i < chain.GetEntries(); i++)
    {
      if (!(i % tenPercent))
      {
        std::cout << i << " / " << chain.GetEntries() << "\t ( " << 10 * i / tenPercent << " %)" << std::endl;
      }
      chain.GetEntry(i);

      outp.cd();
      nt.Fill(e, c);
    }
  }
  outp.cd();
  nt.Write();
  return;
}

int main(int argc, char *argv[])
{
  if (argc != 2)
  {
    std::cout << "Usage: make_trees <event_config_file>" << std::endl;
    return 1;
  }

  std::string configFile(argv[1]);
  std::cout << "Reading from config file " << configFile << std::endl;

  // create the results directory if it doesn't already exist
  std::string outDir;
  ConfigLoader::Open(configFile);
  ConfigLoader::Load("summary", "pruned_ntup_dir", outDir);
  ConfigLoader::Close();

  struct stat st = {0};
  if (stat(outDir.c_str(), &st) == -1)
  {
    mkdir(outDir.c_str(), 0700);
  }

  std::cout << "made " << outDir << std::endl;

  antinufit::EventConfigLoader loader(configFile);
  typedef std::map<std::string, antinufit::EventConfig> EvMap;
  typedef std::vector<std::string> StringVec;
  EvMap toGet = loader.LoadActive();

  // if there is a common path preprend it
  for (EvMap::iterator it = toGet.begin(); it != toGet.end(); ++it)
  {
    const std::string &name = it->first;
    std::cout << "doing " << name << std::endl;
    const std::string &baseDir = it->second.GetNtupBaseDir();
    const StringVec &files = it->second.GetNtupFiles();
    const std::string &outName = it->second.GetPrunedPath();

    std::cout << "Writing from :" << std::endl;
    for (size_t i = 0; i < files.size(); i++)
      std::cout << "\t" << files.at(i) << std::endl;
    std::cout << "to " << outName << std::endl;

    MakeDataSet(files, baseDir, "output", outName);
  }

  return 0;
}
