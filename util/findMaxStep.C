#include <TStyle.h>
#include <TTree.h>

// Script to loop over a bunch of parallel MCMC trees and find which file contained
// the step with the highest LLH. Useful for getting the correct scaled_dist file 
// to input to make_plots

// dataset_dir is level above the fit_dir
void findMaxStep( std::string dataset_dir, std::string fit_dir, int num_files=100 ) {

  // Construct file path
  std::string filepath = std::getenv("DATA_DIR") + dataset_dir + "/" + fit_dir + "/" + fit_dir;

  // Initialise everything to crazy values
  double maxLLH = -999;
  int maxFile = -999;
  int step = -999;
  
  // Loop over files
  for(int i=0; i<num_files; i++){
    
    TString fname = Form("%s_%d/%s_%d_.root",filepath.c_str(),i,fit_dir.c_str(),i);
    std::cout << fname << std::endl;
    
    // Open up the file
    if(!gSystem->AccessPathName(fname)){
      TFile *File = new TFile(fname , "OPEN");

      // Get the chain
      TChain* chain = new TChain("posteriors","");
      chain->Add(fname);

      // Get max LLH
      double llh = chain->GetMaximum("LogL");

      // If LLH > current max LLH
      if(llh > maxLLH){
	// Current file is max file
	maxLLH = llh;
	maxFile = i;
      }
    }
  }

  // Print out max file number
  std::cout << "Max LLH is " << setprecision(10) << maxLLH << " in file " << maxFile << std::endl;

}
