/*
MakeTrees takes in ntuples and prunes out all the branches
we're not interested in so we only have nice lightweight files
to carry around in the fit.
It reads in input files specified in the event config and
literally just loops over events, filling new ntuples
with the quantities we want. These get written to wherever
was specified in the config file.
*/

// Antinu headers
#include <EventConfigLoader.hh>

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

std::vector<std::string> SplitString(std::string str)
{
  std::istringstream buf(str);
  std::istream_iterator<std::string> beg(buf), end;
  std::vector<std::string> tokens(beg, end); // each word of string now in vector
  std::vector<std::string> info;
  std::string dummy = "";

  for (int i = 0; i < tokens.size(); i++)
  { // combine back to reactor name, core number
    if (i == 0)
    {
      dummy = tokens.at(i);
      if (i != tokens.size() - 2)
      {
        dummy += " ";
      }
    }
    else if (i != tokens.size() - 1)
    {
      dummy += tokens.at(i);
      if (i != tokens.size() - 2)
      {
        dummy += " ";
      }
    }
    else
    {
      info.push_back(dummy);
      info.push_back(tokens.at(i));
    }
  }

  return info;
}

void MakeDataSet(const std::vector<std::string> &filenames_,
                 const std::string &baseDir_,
                 const std::string &treeName_,
                 const std::string &outFilename_,
                 std::unordered_map<std::string, int> &reactornameIndex)
{

  // output ntuple
  TFile outp(outFilename_.c_str(), "RECREATE");
  TNtuple nt("pruned", "", "energy:nu_energy:reactorIndex:alphaNClassifier");

  // read the original data
  for (size_t iFile = 0; iFile < filenames_.size(); iFile++)
  {
    std::string fileName = baseDir_ + "/" + filenames_.at(iFile);
    TChain chain(treeName_.c_str());
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

    // read and write
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
        reactorIndex = reactornameIndex[originReactorString];
      }
      else
        reactorIndex = 999;
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
  if (argc != 2)
  {
    std::cout << "Usage: make_trees <event_config_file>" << std::endl;
    return 1;
  }

  std::string configFile(argv[1]);
  std::cout << "Reading from config file " << configFile << std::endl;

  // First read the reactor distance info
  std::unordered_map<std::string, int> reactornameIndex;

  // Read the JSON file
  std::ifstream file("reactors.json");
  if (!file.is_open())
  {
    std::cerr << "Could not open reactors.json file!" << std::endl;
    throw;
  }

  // Parse the JSON file into a json object
  nlohmann::json reactorData;
  file >> reactorData;

  // Loop through the JSON object and fill the maps
  for (const auto &[key, value] : reactorData.items())
  {
    int intKey = std::stoi(key);     // Convert the reactor key to int
    std::string reactorName = value[0]; // First element is the string

    // Fill the map
    reactornameIndex[reactorName] = intKey;
  }

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

    MakeDataSet(files, baseDir, "output", outName, reactornameIndex);
  }

  return 0;
}
