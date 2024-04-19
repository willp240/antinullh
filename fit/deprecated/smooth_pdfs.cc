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
#include <TH1.h>
#include <HistTools.h>
#include <iostream>
using namespace bbfit;

int main(int argc, char *argv[]){
  if (argc != 3){
    std::cout << "\nUsage: .\smooth_pdfs <event_config_file> <pdf_config_file>" << std::endl;
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
	std::string sumDir = pdfDir + "/smoothed_pdfs";
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
		/*
		//smooth psd2	
		std::cout<< "Smoothing PSD2" <<std::endl;
		for(size_t iEBin =0; iEBin< 46; iEBin++){	
			for(size_t iRBin=0; iRBin<6; iRBin++){
				for(size_t iPSD1Bin=0; iPSD1Bin<5; iPSD1Bin++){
					double xx[5] = {0};
					double before = 0;
					for(size_t iPSD2Bin=0; iPSD2Bin<5; iPSD2Bin++){
						std::vector<size_t> indexVec;
						indexVec.push_back(iEBin);
						indexVec.push_back(iRBin);
						indexVec.push_back(iPSD1Bin);
						indexVec.push_back(iPSD2Bin);
						size_t index = dist.FlattenIndices(indexVec);
						xx[iPSD2Bin] = dist.GetBinContent(index);
						before+=xx[iPSD2Bin];
					}
					if (!before) continue;
					TH1::SmoothArray(5, xx, 1);
				
					//replace
					double after = 0;
					for(size_t iPSD2Bin=0; iPSD2Bin<5; iPSD2Bin++){
						std::vector<size_t> indexVec;
						indexVec.push_back(iEBin);
						indexVec.push_back(iRBin);
						indexVec.push_back(iPSD1Bin);
						indexVec.push_back(iPSD2Bin);
						size_t index = dist.FlattenIndices(indexVec);
						dist.SetBinContent(index, xx[iPSD2Bin]);
						after+=xx[iPSD2Bin];
					}
					if ((!before) || (!after)) continue;
					dist.Scale(before/after);
				}
			}
		}
		std::cout<< "Integral after: " << dist.Integral()<<std::endl;


		
		//smooth psd1
		std::cout<< "Smoothing PSD1" <<std::endl;
		for(size_t iEBin =0; iEBin< 46; iEBin++){	
			for(size_t iRBin=0; iRBin<6; iRBin++){
				for(size_t iPSD2Bin=0; iPSD2Bin<5; iPSD2Bin++){
					double xx[5] = {0};
					double before = 0;
					for(size_t iPSD1Bin=0; iPSD1Bin<5; iPSD1Bin++){
						std::vector<size_t> indexVec;
						indexVec.push_back(iEBin);
						indexVec.push_back(iRBin);
						indexVec.push_back(iPSD1Bin);
						indexVec.push_back(iPSD2Bin);
						size_t index = dist.FlattenIndices(indexVec);
						xx[iPSD1Bin] = dist.GetBinContent(index);
						before+=xx[iPSD1Bin];
					}
					if (!before) continue;
					TH1::SmoothArray(5, xx, 1);
				
					//replace
					double after = 0;
					for(size_t iPSD1Bin=0; iPSD1Bin<5; iPSD1Bin++){
						std::vector<size_t> indexVec;
						indexVec.push_back(iEBin);
						indexVec.push_back(iRBin);
						indexVec.push_back(iPSD1Bin);
						indexVec.push_back(iPSD2Bin);
						size_t index = dist.FlattenIndices(indexVec);
						dist.SetBinContent(index, xx[iPSD1Bin]);
						after+=xx[iPSD1Bin];
					}
					if ((!before) || (!after)) continue;
					dist.Scale(before/after);
				}
			}
		}
		std::cout<< "Integral after: " << dist.Integral()<<std::endl;
		*/
		//smooth E	
		std::cout<< "Smoothing E" <<std::endl;
		for(size_t iRBin=0; iRBin<6; iRBin++){
			for(size_t iPSD1Bin=0; iPSD1Bin<5; iPSD1Bin++){
				for(size_t iPSD2Bin=0; iPSD2Bin<5; iPSD2Bin++){
					double xx[48] = {0};
					double before = 0;
					for(size_t iEBin =0; iEBin< 48; iEBin++){	
						std::vector<size_t> indexVec;
						indexVec.push_back(iEBin);
						indexVec.push_back(iRBin);
						indexVec.push_back(iPSD1Bin);
						indexVec.push_back(iPSD2Bin);
						//std::cout<< "Vec: " << indexVec.size() << std::endl;
						size_t index = dist.FlattenIndices(indexVec);
						//std::cout<< "Index: " << index<< std::endl;
						//std::cout<< "Element before: " << dist.GetBinContent(3)<< std::endl;
						//std::cout<< "Element before: " << dist.GetBinContent(index)<< std::endl;
						//std::cout<< "Element before: " << xx[iEBin]<< std::endl;
						xx[iEBin] = dist.GetBinContent(index);
						//std::cout<< "Element before: " << xx[iEBin]<< std::endl;
						before+=xx[iEBin];
					}
					if (!before) continue;
					
					TH1::SmoothArray(48, xx, 1);
				
					//replace
					double after = 0;
					for(size_t iEBin =0; iEBin< 48; iEBin++){	
						std::vector<size_t> indexVec;
						indexVec.push_back(iEBin);
						indexVec.push_back(iRBin);
						indexVec.push_back(iPSD1Bin);
						indexVec.push_back(iPSD2Bin);
						size_t index = dist.FlattenIndices(indexVec);
						
						//std::cout<< "EBin " << iEBin << std::endl;
						//std::cout<< "Element after " << dist.GetBinContent(index) << std::endl;
						// //std::cout<< "Element after " << xx[iEBin] << std::endl;
						dist.SetBinContent(index, xx[iEBin]);
						//std::cout<< "Element after " << dist.GetBinContent(index) << std::endl;
						after+=xx[iEBin];
					}
					if ((!before) || (!after)) continue;
					dist.Scale(before/after);
				}
			}
		}
		std::cout<< "Integral after: " << dist.Integral()<<std::endl;


		
		//smooth R
/*		std::cout<< "Smoothing R" <<std::endl;
		
		for(size_t iPSD1Bin=0; iPSD1Bin<5; iPSD1Bin++){
			for(size_t iPSD2Bin=0; iPSD2Bin<5; iPSD2Bin++){
				for(size_t iEBin =0; iEBin< 46; iEBin++){	
					double xx[6] = {0};
					double before = 0;
					for(size_t iRBin=0; iRBin<6; iRBin++){
						std::vector<size_t> indexVec;
						indexVec.push_back(iEBin);
						indexVec.push_back(iRBin);
						indexVec.push_back(iPSD1Bin);
						indexVec.push_back(iPSD2Bin);
						size_t index = dist.FlattenIndices(indexVec);
						xx[iRBin] = dist.GetBinContent(index);
						before+=xx[iRBin];
					}
					if (!before) continue;

					TH1::SmoothArray(6, xx, 1);
				
					//replace
					double after = 0;
					for(size_t iRBin =0; iRBin< 6; iRBin++){	
						std::vector<size_t> indexVec;
						indexVec.push_back(iEBin);
						indexVec.push_back(iRBin);
						indexVec.push_back(iPSD1Bin);
						indexVec.push_back(iPSD2Bin);
						size_t index = dist.FlattenIndices(indexVec);
						dist.SetBinContent(index, xx[iRBin]);
						after+=xx[iRBin];
					}
					if ((!before) || (!after)) continue;
					dist.Scale(before/after);
				}
			}
		}
		std::cout<< "Integral after: " << dist.Integral()<<std::endl;
*/	
	
		
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
