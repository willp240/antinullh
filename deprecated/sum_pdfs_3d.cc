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
#include <AxisCollection.h>
#include <BinAxis.h>

#include <HistTools.h>
#include <iostream>
using namespace bbfit;

int main(int argc, char *argv[]){
  if (argc != 3){
    std::cout << "\nUsage: .\sum_pdfs <event_config_file> <pdf_config_file>" << std::endl;
    return 1;
  }
    
  std::string evConfigFile(argv[1]);
  std::string pdfConfigFile(argv[2]);


  std::cout << "\nReading from config files: "   << std::endl
	    << "\t" << evConfigFile << ",\n "  
	    << "\t" << pdfConfigFile
	    << std::endl;

    
  // load up the pdfs
  DistConfigLoader pLoader(pdfConfigFile);
  DistConfig pConfig = pLoader.Load();
  std::string pdfDir = pConfig.GetPDFDir();
	std::string sumDir = pdfDir + "/summed_composite_pdfs_full_range";
  struct stat st = {0};
  if (stat(sumDir.c_str(), &st) == -1) {
    mkdir(sumDir.c_str(), 0700);
  }

  std::cout << "\nSaving pdfs to " << sumDir << std::endl;

  // and another one for the projections - there will be loads
  std::string projDir = sumDir + "/projections";
  if (stat(projDir.c_str(), &st) == -1) {
    mkdir(projDir.c_str(), 0700);
  }
  
  std::cout << "\nSaving projections logs to " << projDir << std::endl;

  // load up all the event types we want pdfs for
	std::vector<BinnedED> dists;
	
	
  typedef std::map<std::string, EventConfig> EvMap;
  EventConfigLoader loader(evConfigFile);
  EvMap toGet = loader.LoadActive();

 
  // now make and fill the pdfs
  for(EvMap::iterator it = toGet.begin(); it != toGet.end(); ++it){    
    std::cout << "Retrieving distribution for " << it->first << std::endl;
    std::string distPath = pdfDir + "/" + it->first + ".h5";
		
    BinnedED dist = BinnedED(it->first, IO::LoadHistogram(distPath));
    dists.push_back(dist);
    
    
    		
    //full range
    double energyLowPSD = 1.8;
    double energyHighPSD = 3.0; 
    

    //in ROI, ignore running of PSD with energy in radial slices 
    //get the boundary energy bins
    AxisCollection axCol = dist.GetAxes();
    BinAxis axE = axCol.GetAxis(0);
    size_t lowEBin = axE.FindBin(energyLowPSD);
    size_t highEBin = axE.FindBin(energyHighPSD);
    
    std::vector<std::map<int, double> > vectorOfMaps;
    
    for(size_t iRBin=0; iRBin<6; iRBin++){
      std::map<int, double> contentMap;
		  
      double integralPSD = 0;
      for(size_t iPSD1Bin=0; iPSD1Bin<5; iPSD1Bin++){
	
	
	//get contents and summ
	double summedBinContent=0;
	double summedBins =0;
	//only in roi
	for(size_t iEBin =0; iEBin< 48; iEBin++){
	  //if ((iEBin<=lowEBin) || (iEBin>=highEBin)) continue;
	  
	  std::vector<size_t> indexVec;
	  indexVec.push_back(iEBin);
	  indexVec.push_back(iRBin);
	  indexVec.push_back(iPSD1Bin);
	  
	  size_t index = dist.FlattenIndices(indexVec);
	  summedBinContent += dist.GetBinContent(index);
	  summedBins++;
	}
	summedBinContent/=summedBins;
	integralPSD+=summedBinContent;
	std::cout<< "summedBinContent: " << summedBinContent << std::endl;
	std::cout<< "summedBins: " << summedBins <<std::endl;
				
	contentMap.insert(std::make_pair(iPSD1Bin, summedBinContent));
	
      }
      
      std::map<int, double>::iterator it = contentMap.begin();
      std::map<int, double> contentMapNormalised;
      while (it!=contentMap.end()){
	std::cout<< "it->first: " << it->first << "  it->second: " << it->second <<std::endl;
	if (integralPSD == 0){
	  contentMapNormalised.insert(std::make_pair(it->first, 0.04));	
	  std::cout<< "it->first: " << it->first << "  it->second norm: " << 0.04 <<std::endl;
	}else{
	  contentMapNormalised.insert(std::make_pair(it->first, (it->second)/integralPSD));	
	  std::cout<< "it->first: " << it->first << "  it->second norm: " << (it->second)/integralPSD <<std::endl;
	}
	it++;
      }
			
      vectorOfMaps.push_back(contentMapNormalised);
      
    }
    
		
    //sum across PSD in energy outside of ROI
    //get the boundary energy bins
    /*AxisCollection axCol = dist.GetAxes();
      BinAxis axE = axCol.GetAxis(0);
      size_t lowEBin = axE.FindBin(energyLowPSD);
		size_t highEBin = axE.FindBin(energyHighPSD);
    */
    
    for(size_t iEBin =0; iEBin< 48; iEBin++){
      //if ((iEBin>lowEBin) && (iEBin<highEBin)) continue;
      
      for(size_t iRBin=0; iRBin<6; iRBin++){

	
	//get contents and summ
	double summedBinContent=0;
	for(size_t iPSD1Bin=0; iPSD1Bin<5; iPSD1Bin++){
	  
	  std::vector<size_t> indexVec;
	  indexVec.push_back(iEBin);
	  indexVec.push_back(iRBin);
	  indexVec.push_back(iPSD1Bin);
						
	  size_t index = dist.FlattenIndices(indexVec);
	  summedBinContent += dist.GetBinContent(index);
	  
	}
	summedBinContent /=5;
	
	//replace
	for(size_t iPSD1Bin=0; iPSD1Bin<5; iPSD1Bin++){
	  
	  std::vector<size_t> indexVec;
	  indexVec.push_back(iEBin);
	  indexVec.push_back(iRBin);
	  indexVec.push_back(iPSD1Bin);
	  
	  size_t index = dist.FlattenIndices(indexVec);
	  dist.SetBinContent(index, summedBinContent);
	}
	
	
      }
    }
    
    // normalise 
    std::cout<< "Integral ER: " << dist.Integral() << std::endl;	
    if(dist.Integral()){
      std::cout<< "Normalising" << std::endl;	
      dist.Normalise();
    }
    
    //replace
    for(size_t iRBin=0; iRBin<6; iRBin++){
      std::map<int, double> contentMap = vectorOfMaps.at(iRBin);
      double normCheck = 0.0;
      for(size_t iPSD1Bin=0; iPSD1Bin<5; iPSD1Bin++){
	
	double toReplace = contentMap.find(iPSD1Bin)->second;
	std::cout<< "toReplace: " << toReplace<< std::endl;
	normCheck+=toReplace;
					
	for(size_t iEBin =0; iEBin< 48; iEBin++){
	  //if ((iEBin<=lowEBin) && (iEBin>=highEBin)) continue;
	  std::vector<size_t> indexVec;
	  indexVec.push_back(iEBin);
	  indexVec.push_back(iRBin);
	  indexVec.push_back(iPSD1Bin);
	  
	  size_t index = dist.FlattenIndices(indexVec);
	  double oldContent = dist.GetBinContent(index);
	  dist.SetBinContent(index, (oldContent*toReplace));
	}
	
	
      }
      std::cout<< "NormCheck: "<< normCheck <<std::endl;
    }
    
    
    // normalise 
    std::cout<< "Integral" << dist.Integral() << std::endl;	
    if(dist.Integral()){
      std::cout<< "Normalising" << std::endl;	
      dist.Normalise();
    }
	
    
    // detect zero bins
    int zeroBins = 0;
    for(int i = 0; i<dist.GetNBins(); i++){
			if (!dist.GetBinContent(i)){
			  //std::cout << "content "<< dist.GetBinContent(i) << " for bin " << i << std::endl;
			  zeroBins++;
			}
    }
    std::cout<< "Zero bins in dist " << it->first << " is " << zeroBins << std::endl;
    
    // save as h5
    IO::SaveHistogram(dist.GetHistogram(), sumDir + "/" + it->first + ".h5");
    
    // save as a root histogram if possible
    if(dist.GetNDims() <= 2)
        IO::SaveHistogram(dist.GetHistogram(), sumDir + "/" + it->first + ".root");
    
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
