
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
const double rav = 6005;

double 
Reff(double x_, double y_, double z_){
  return pow((sqrt(x_ * x_ + y_ * y_ + z_ * z_)/rav), 3);
}

void
MakeDataSet(const std::vector<std::string>& filenames_, 
	    const std::string& baseDir_, 
	    const std::string& treeName_,
	    const std::string& outFilename_
	    ){

  // output ntuple
  TFile outp(outFilename_.c_str(), "RECREATE");
  TNtuple nt("pruned", "", "energy:fitValid:reff:qmcdep:bipoCumul:biPoLikelihood214:itr:timePSD:anglePSD");
 
  // read the original data
  for(size_t iFile = 0; iFile < filenames_.size(); iFile++){
    std::string fileName = baseDir_ + "/" + filenames_.at(iFile);
    TChain c(treeName_.c_str());
    c.Add(fileName.c_str());
    std::cout << fileName << "\t" << c.GetEntries() << "  entries"<< std::endl;
    
    int tenPercent = c.GetEntries()/10;

    Double_t e;
    Bool_t   v;
    Double_t x;
    Double_t y;
    Double_t z;
    Double_t mce;
    Double_t bpL214;
    Double_t bpCumul;
    Double_t itr;
		Double_t timePSD;
		Double_t anglePSD;

    c.SetBranchAddress("energy", &e);
    c.SetBranchAddress("fitValid", &v);
    c.SetBranchAddress("posx", &x);
    c.SetBranchAddress("posy", &y);
    c.SetBranchAddress("posz", &z);
    c.SetBranchAddress("mcEdepQuenched", &mce);
    c.SetBranchAddress("biPoCumul", &bpCumul);
    c.SetBranchAddress("biPoLikelihood214", &bpL214);
    c.SetBranchAddress("itr", &itr);
		c.SetBranchAddress("ext0NuTimeTl208AVNaive", &timePSD);
		c.SetBranchAddress("ext0NuAngleTl208AV", &anglePSD);
    
    // read and write
    Double_t reff;
    
    for(int i = 0; i < c.GetEntries(); i++){
      if(!(i % tenPercent)){
	std::cout << i << " / " << c.GetEntries() << "\t ( " << 10 * i/tenPercent << " %)" << std::endl;
      }
      c.GetEntry(i);
      reff = Reff(x, y, z);
      
      outp.cd();
      nt.Fill(e, v, reff, mce, bpCumul, bpL214, itr, timePSD, anglePSD);
    }
  }
  outp.cd();
  nt.Write();
  return;
}

int main(int argc, char *argv[]){
  if (argc != 2){
    std::cout << "Usage: make_trees <event_config_file>" << std::endl;
    return 1;
  }
    
  std::string configFile(argv[1]);
  std::cout << "Reading from config file "  << configFile << std::endl;

  // create the results directory if it doesn't already exist
  std::string outDir;
  ConfigLoader::Open(configFile);
  ConfigLoader::Load("summary", "pruned_ntup_dir", outDir);
  ConfigLoader::Close();
  
  struct stat st = {0};
  if (stat(outDir.c_str(), &st) == -1) {
    mkdir(outDir.c_str(), 0700);
  }

  bbfit::EventConfigLoader loader(configFile);
  typedef std::map<std::string, bbfit::EventConfig> EvMap;
  typedef std::vector<std::string> StringVec;
  EvMap toGet = loader.LoadActive();

  // if there is a common path preprend it
  for(EvMap::iterator it = toGet.begin(); it != toGet.end(); ++it){
    const std::string& name = it->first;
    const std::string& baseDir = it->second.GetNtupBaseDir();  
    const StringVec& files = it->second.GetNtupFiles();
    const std::string& outName = it->second.GetPrunedPath();

    std::cout << "Writing from :" << std::endl;
    for(size_t i = 0; i < files.size(); i++)
      std::cout << "\t" << files.at(i) << std::endl;
    std::cout << "to " << outName << std::endl;

    MakeDataSet(files, baseDir, "output", outName);
  }
    
  return 0;
}
