#include <DistConfigLoader.hh>
#include <DistConfig.hh>
#include <DistBuilder.hh>
#include <EventConfig.hh>
#include <EventConfigLoader.hh>
#include <iostream>
#include <sys/stat.h>
#include <BinnedED.h>
#include <DistFiller.h>
#include <ROOTNtuple.h>
#include <IO.h>
#include <DistTools.h>
#include <TH1D.h>
#include <TH2D.h>
#include <HistTools.h>
#include <iostream>
using namespace antinufit;

int main(int argc, char *argv[]){
  if (argc != 4){
    std::cout << "\nUsage: make_trees <event_config_file> <pdf_config_file>" << std::endl;
    return 1;
  }
    
  std::string evConfigFile(argv[1]);
  std::string pdfConfigFile(argv[2]);


  std::cout << "\nReading from config files: "   << std::endl
	    << "\t" << evConfigFile << ",\n "  
	    << "\t" << pdfConfigFile << ",\n"
	    << std::endl;

    
  // load up the pdf configuration data too
  DistConfigLoader pLoader(pdfConfigFile);
  DistConfig pConfig = pLoader.Load();
 
  // create the pdf directory if it doesn't already exist
  std::string pdfDir = pConfig.GetPDFDir();
  struct stat st = {0};
  if (stat(pdfDir.c_str(), &st) == -1) {
    mkdir(pdfDir.c_str(), 0700);
  }

  std::cout << "\nSaving pdf logs to " << pdfDir << std::endl;

  // and another one for the projections - there will be loads
  std::string projDir = pdfDir + "/projections";
  if (stat(projDir.c_str(), &st) == -1) {
    mkdir(projDir.c_str(), 0700);
  }
  
  std::cout << "\nSaving projections logs to " << projDir << std::endl;

  // load up all the event types we want pdfs for
  typedef std::map<std::string, EventConfig> EvMap;
  EventConfigLoader loader(evConfigFile);
  EvMap toGet = loader.LoadActive();
 
  // now make and fill the pdfs
  for(EvMap::iterator it = toGet.begin(); it != toGet.end(); ++it){    
    std::cout << "Building distribution for " << it->first << std::endl;
    
    // find the dataset
    DataSet* dataSet;
    try{
        dataSet = new ROOTNtuple(it->second.GetSplitPdfPath(), "pruned");
    }
    catch(const IOError& e_){
        std::cout << "Warning: skipping " << it-> first << " couldn't open data set:\n\t" << e_.what() << std::endl;
        continue;
    }

    // create and fill
    BinnedED dist = DistBuilder::Build(it->first, pConfig, dataSet, log);

    std::vector<double> gen;
    gen.push_back(dist.Integral());

    delete dataSet;

    // normalise
    if(dist.Integral())
      dist.Normalise();

    // save a copy of the cut log
    log.SaveAs(it->first, pdfDir + "/" + it->first + ".txt");

    // save as h5
    IO::SaveHistogram(dist.GetHistogram(), pdfDir + "/" + it->first + ".h5");

    // save as a root histogram if possible
    if(dist.GetNDims() <= 2)
        IO::SaveHistogram(dist.GetHistogram(), pdfDir + "/" + it->first + ".root");

    TFile *fOut = TFile::Open((pdfDir + "/" + it->first + ".root").c_str(), "UPDATE");
    fOut->WriteObject(&gen, "nGeneratedEvents");
    fOut->Close();

    // HigherD save the projections
    if(dist.GetNDims() > 1){
      std::vector<BinnedED> projs = HistTools::GetVisualisableProjections(dist);

      // save them as apropriate
      for(size_t i = 0; i < projs.size(); i++){
          const BinnedED& proj = projs.at(i);

          if(projs.at(i).GetNDims() == 1)
              DistTools::ToTH1D(proj).SaveAs((projDir + "/" + proj.GetName() + ".root").c_str());
          else
              DistTools::ToTH2D(proj).SaveAs((projDir + "/" + proj.GetName() + ".root").c_str());
      }
    }
  }

  return 0;
}
