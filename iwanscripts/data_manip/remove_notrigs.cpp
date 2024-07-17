#include <TTree.h>
#include <TChain.h>
#include <TVector3.h>
#include <RAT/DS/UniversalTime.hh>
#include <iostream>
#include <list>
#include <cmath>
#include <TFile.h>
#include <TH1D.h>
#include <sstream>

using RAT::DS::UniversalTime;

void
Clean(const std::string& infiles_, const std::string& outfile_, int testing_){
    TChain c("output");
    c.Add(infiles_.c_str());

    TFile output(outfile_.c_str(), "RECREATE");
    TTree* newTree = c.CloneTree(0);

    Int_t evIndex;
    c.SetBranchAddress("evIndex", &evIndex);

    int onePercent = c.GetEntries()/100;
    for(int i = 0; i < c.GetEntries(); i++){
      c.GetEntry(i);
      if(!(i%onePercent))
	std::cout << i/onePercent << "% done " << std::endl;

      bool passes = evIndex > -1;

      if(testing_ && i > 10)
	break;

      // save that bad boy
      if(passes)
	newTree->Fill();
    }
    output.cd();
    newTree->Write();
    output.Close();
}    



int main(int argc, char* argv[]){
  if(argc != 4){
    std::cout << "Usage: ./remove_notrigs <infiles> <outfile> <testing>" << std::endl;
    return 1;
  }

  int  testing;
  std::istringstream(argv[3]) >> testing;

  Clean(std::string(argv[1]), std::string(argv[2]), testing);
  return 0;
}
