#include <DistConfigLoader.hh>
#include <DistConfig.hh>
#include <DistBuilder.hh>
#include <iostream>
#include <sys/stat.h>
#include <BinnedED.h>
#include <DistFiller.h>
#include <ROOTNtuple.h>
#include <IO.h>
#include <DistTools.h>
#include <TH1D.h>
#include <TH2D.h>
#include <FitConfig.hh>
#include <FitConfigLoader.hh>
#include <HistTools.h>
#include <iostream>
using namespace bbfit;

int main(int argc, char *argv[]){
  if (argc != 3){
    std::cout << "\nUsage: slice_pdfs <fit_config_file> <pdf_config_file>" << std::endl;
    return 1;
  }
    
  std::string fitConfigFile(argv[1]);
  std::string pdfConfigFile(argv[2]);


  std::cout << "\nReading from config files: "   << std::endl
	    << "\t" << fitConfigFile << ",\n "  
	    << "\t" << pdfConfigFile 
	    << std::endl;

    
  // load up the pdf configuration data too
  DistConfigLoader pLoader(pdfConfigFile);
  DistConfig pConfig = pLoader.Load();
  std::string pdfDir = pConfig.GetPDFDir();
	
  // create directory for the slices - there will be loads
	struct stat st = {0};
  std::string sliceDir = pdfDir + "/slices";
  if (stat(sliceDir.c_str(), &st) == -1) {
    mkdir(sliceDir.c_str(), 0700);
  }
  
  std::cout << "\nSaving slices to " << sliceDir << std::endl;


  // want slices for the pdfs that are use in the fit, i.e. in mcmcconfig
	FitConfigLoader mcLoader(fitConfigFile);
	FitConfig fConfig = mcLoader.LoadActive();
  typedef std::set<std::string> StringSet;
  StringSet distsToFit = fConfig.GetParamNames();

	
  for(StringSet::iterator it = distsToFit.begin(); it != distsToFit.end();
      ++it){
    std::string distPath = pdfDir + "/" + *it + ".h5";
    BinnedED dist = BinnedED(*it, IO::LoadHistogram(distPath));
		dist.SetObservables(pConfig.GetBranchNames());
		
    // only get slices for higher D
    if(dist.GetNDims() > 1){
      std::vector<BinnedED> slices = HistTools::Get1DSlices(dist, "energy");

      // save them as apropriate
      for(size_t i = 0; i < slices.size(); i++){
          const BinnedED& slice = slices.at(i);
					DistTools::ToTH1D(slice).SaveAs((sliceDir + "/" + slice.GetName() + ".root").c_str());
          
      }
    }
  }

  return 0;
}
