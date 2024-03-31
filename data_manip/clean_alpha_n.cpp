#include <TTree.h>
#include <TChain.h>
#include <TVector3.h>
#include <RAT/DS/UniversalTime.hh>
#include <iostream>
#include <cmath>
#include <TFile.h>
#include <TH1D.h>
#include <sstream>
#include <sys/stat.h>

using RAT::DS::UniversalTime;

void
Clean(const std::string& infiles_, const std::string& outfile_, double eMin_, double deltaT_, const std::string histDir_, int testing){
    TChain c("output");
    c.Add(infiles_.c_str());

    TFile output(outfile_.c_str(), "RECREATE");
    TTree* newTree = c.CloneTree(0);
    
    struct stat st = {0};
    if (stat(histDir_.c_str(), &st) == -1) {
      mkdir(histDir_.c_str(), 0700);
    }
    std::cout << "in here " << std::endl;
    

    Int_t secs;
    Int_t nsecs;
    Int_t days;
    Int_t evIndex;

    Double_t energy;

    Double_t x;
    Double_t y;
    Double_t z;
    
    Bool_t fitValid;
    Bool_t nextFitValid;

    Int_t nhit;
    c.SetBranchAddress("uTDays", &days);
    c.SetBranchAddress("uTSecs", &secs);
    c.SetBranchAddress("uTNSecs", &nsecs);
    c.SetBranchAddress("energy", &energy);
    c.SetBranchAddress("fitValid", &fitValid);
    c.SetBranchAddress("evIndex", &evIndex);
    c.SetBranchAddress("posx", &x);
    c.SetBranchAddress("posy", &y);
    c.SetBranchAddress("posz", &z);
    c.SetBranchAddress("nhits", &nhit);

    Int_t nextDays;
    Int_t nextSec;
    Int_t nextNSec;

    Double_t nextEnergy;
    
    int roiCount = 0;
    int roiPassCount = 0;
    
    bool flag1 = true; // assume that the first event is invalidated by an event we didn't catch
    bool flag2;

    TH1D beforeEnergy("b", "", 1000, 0, 10);
    TH1D afterEnergy("a", "", 1000, 0, 10);

    
    int onePercent = c.GetEntries()/100;
    for(int i = 0; i < c.GetEntries(); i++){
      if(onePercent && !(i%onePercent))
	std::cout << i/onePercent << "% done " << std::endl;

      if(testing && i > 10)
	break;

      // assume the last event is invalidated by one after we didn't catch
      if(i == (c.GetEntries() -1)){
	flag2 = true;
      }
      else{
	// get the next event
	c.GetEntry(i+1);
	nextDays = days;
	nextSec = secs;
	nextNSec = nsecs;
	nextEnergy = energy;
	nextFitValid = fitValid;

	// and this one
	c.GetEntry(i);        
	
	bool isCoincidence;
	if(std::abs(nextDays - days) > 0)
	  isCoincidence = false;
	
	else if(std::abs(nextSec - secs) > 0)
	  isCoincidence = false;
	
	else if(std::abs(nextNSec - nsecs) < deltaT_)
	  isCoincidence = true;
	
	else
	  isCoincidence = false;     
	
	flag2 = (isCoincidence && (energy > eMin_) && (nextEnergy > eMin_) && fitValid && nextFitValid);
      }
      
      bool passes = !(flag1 || flag2);

      beforeEnergy.Fill(energy);
      if(passes)
	afterEnergy.Fill(energy);

      if((energy > 2.47) && (energy < 2.64) && (TVector3(x, y, z).Mag() < 3500) && fitValid == 1){
	roiCount++;
	if(passes)
	  roiPassCount++;
      }
      
      
      // save that bad boy
      if(passes)
	newTree->Fill();
      
      // shuffle the flag along, if this event is tagged by it's coincidence with the next (flag2)
      // the next one is flagged by its coincidence with the one before (flag1)
      flag1 = flag2;
    }

    std::cout << "For the roi  " << roiPassCount << " pass out of " << roiCount << std::endl;
    std::cout << c.GetEntries() << std::endl;

    beforeEnergy.SaveAs((histDir_ + "/an_energy_before.root").c_str());
    afterEnergy.SaveAs((histDir_ + "/an_energy_after.root").c_str());

    newTree->Write();
    output.Close();
}    



int main(int argc, char* argv[]){
  if(argc != 7){
    std::cout << "Usage: ./clean_alpha_n <infiles> <outfile> e_min time_diff hist_dir testing" << std::endl;
    return 1;
  }
    
  float energy;
  std::istringstream(argv[3]) >> energy;

  float timeDiff;
  std::istringstream(argv[4]) >> timeDiff;

  int testing;
  std::istringstream(argv[6]) >> testing;

  Clean(std::string(argv[1]), std::string(argv[2]), energy, timeDiff, std::string(argv[5]), testing);
  return 0;
}
