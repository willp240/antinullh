#include <EventConfigLoader.hh>
#include <EventConfig.hh>
#include <iostream>
#include <IO.h>
#include <ROOTNtuple.h>
#include <DataSetGenerator.h>
#include <OXSXDataSet.h>
#include <sstream>
#include <fstream>
#include <Formatter.hpp>
#include <ConfigLoader.hh>
#include <sys/stat.h>
#include <Rand.h>
using namespace bbfit;

void SaveRemainders(DataSetGenerator &dsGen, std::vector<DataSet*> &dataSets, std::vector<std::string> &names, const std::string &configFile_){

  std::string outDirPdf;
  ConfigLoader::Open(configFile_);
  ConfigLoader::Load("summary", "split_ntup_dir_pdf", outDirPdf);
  ConfigLoader::Close();

  struct stat st = {0};
	if (stat(outDirPdf.c_str(), &st) == -1) {
    mkdir(outDirPdf.c_str(), 0700);
  }
	
  // save what's left over as independent data sets
  OXSXDataSet* remainder = NULL;
  int countsTaken = 0;
  std::vector<int> content;
  for(int iSet = 0; iSet < dataSets.size(); iSet++){
    std::cout << "Assembling the remainder for  " << names.at(iSet) << std::endl;

    remainder = dsGen.AllRemainingEvents(iSet, &countsTaken);
    content.push_back(countsTaken);
    
    std::cout << "\t.. and saving" << std::endl;
    IO::SaveDataSet(*remainder, Formatter() << outDirPdf << "/" << names.at(iSet) << ".root", "pruned");		
    delete remainder;
  }

  std::ofstream fs;
  fs.open((outDirPdf + "/mc_remainders.txt").c_str());
  for(size_t i = 0; i < names.size(); i++)
    fs << names.at(i) << "\t" << content.at(i) << "\n";
  fs.close();  
  
  std::cout << "\n\n Left over events written to  " << outDirPdf
	    << "\t with logfile " << outDirPdf + "mc_remainders.txt" << std::endl;    
}

void
MakeDataSets(const std::string& configFile_, double liveTime_, int nDataSets_, bool replaceEvents_){
  // load up the rates of the different event types
  EventConfigLoader loader(configFile_);

  typedef std::map<std::string, EventConfig>  EvMap;
  EvMap active = loader.LoadActive();
 
  // pull out the useful parts
  std::vector<std::string> names;
  std::vector<double> rates;
  std::vector<DataSet*> dataSets;
  std::vector<bool> flags;
  std::vector<bool> bootstraps;
	
  for(EvMap::iterator it = active.begin(); it != active.end(); ++it){
    DataSet* ds = new ROOTNtuple(it->second.GetPrunedPath(), "pruned");
    // note the correction for the number of events generated e.g scintEdep cut
    double expectedCounts = it->second.GetRate() * liveTime_;
    
    if(it->second.GetNGenerated())
      expectedCounts *= ds->GetNEntries()/double(it->second.GetNGenerated());

    if(!expectedCounts)
      std::cout << "\n(" << ds->GetNEntries() << "\t" << it->second.GetNGenerated() << "\t" << liveTime_ << "\t" << it->second.GetRate() << ")\n" << std::endl;

    if(!ds->GetNEntries()){
      std::cout << "Warning:: skipping " << it->first << "  no events to choose from" << std::endl;
      continue;
    }

    //Rand::SetSeed(0);

    dataSets.push_back(ds);
    names.push_back(it->first);
    rates.push_back(expectedCounts);
    flags.push_back(!it->second.GetRandomSplit());
    bootstraps.push_back(replaceEvents_);

    //Could use the below instead of the line above to make sure to just replace events with low stats
    //This creates somewhat semi-independent datasets though!!!
    /*
    std::cout<< "Name: " << (it->first).c_str() <<std::endl;
    if ( ((it->first)=="tl208_av") ||
	 ((it->first)=="tl208_hup") ||
				 ((it->first)=="tl208_hdr") ||
				 ((it->first)=="tl208_exwater") ||
				 ((it->first)=="bi214_av") ||
				 ((it->first)=="bi214_od") ||
				 ((it->first)=="bi214_exwater")
			){
			bootstraps.push_back(1);
			std::cout<< "Overring bootstrap! " <<std::endl;
		
		}else{
			bootstraps.push_back(replaceEvents_);
		}

		std::cout<< "Name: " << it->first <<std::endl;
    */

    std::cout << "Loading data set " 
	      << it->second.GetPrunedPath() 
	      << " with expected counts "
	      << expectedCounts
	      << "\n";

	}

  DataSetGenerator dsGen;
  dsGen.SetDataSets(dataSets);
  dsGen.SetExpectedRates(rates);
  dsGen.SetSequentialFlags(flags);
  dsGen.SetBootstrap(bootstraps);
		
	
  std::string outDirFake;
  ConfigLoader::Open(configFile_);
  ConfigLoader::Load("summary", "split_ntup_dir_fake", outDirFake);
  ConfigLoader::Close();

  struct stat st = {0};
  if (stat(outDirFake.c_str(), &st) == -1) {
    mkdir(outDirFake.c_str(), 0700);
  }

  std::cout << "Generating " << nDataSets_ << " data sets with livetime " << liveTime_
	    << "  including poisson fluctuations...\n" << std::endl;

  // actually generate the events
  std::vector<int> content;
  for(int iSet = 0; iSet < nDataSets_; iSet++){
    std::cout << "DataSet #" << iSet << std::endl;
    std::string outPath = Formatter() << outDirFake << "/fake_data_lt_" << liveTime_ << "__" << iSet;

    //OXSXDataSet ds = dsGen.ExpectedRatesDataSet(&content);
    OXSXDataSet ds = dsGen.PoissonFluctuatedDataSet(&content);
    std::ofstream fs;
    fs.open((outPath + ".txt").c_str());
    for(size_t i = 0; i < names.size(); i++)
      fs << names.at(i) << "\t" << content.at(i) << "\n";
    fs.close();

    // save to a new ROOT tree
    IO::SaveDataSet(ds, outPath + ".root", "pruned");

    std::cout << "\t .. written " << ds.GetNEntries() << " events to "  << outPath + ".root"
	      << "\t with logfile " << outPath + ".txt\n\n" << std::endl;
  }

	//If events not replaced, save the remainders, otherwise it doesn't make sense
	if (!replaceEvents_)
		SaveRemainders(dsGen, dataSets, names, configFile_);

}


int main(int argc, char *argv[]){
  if(argc != 5){
    std::cout << "\nUsage: ./split_data <event_config_file> <livetime> <n_data_sets> <replace_events(0 or 1)>" << std::endl;
    return 1;
  }

  std::string configFile(argv[1]);
  int nDataSets;
  std::istringstream(argv[3]) >> nDataSets;

  double liveTime;
  std::istringstream(argv[2]) >> liveTime;

	bool replaceEvents;
	std::istringstream(argv[4]) >> replaceEvents;

	MakeDataSets(configFile, liveTime, nDataSets, replaceEvents);
   
  return 0;
}
