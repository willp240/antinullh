#include <BayesIntervalCalc.h>
#include <DistTools.h>
#include <TFile.h>
#include <TH1D.h>
#include <IO.h>
#include <iostream>
#include <string>
#include <sstream>

int main(int argc, char *argv[]){
    if(argc != 4){
        std::cout << "Usage:\n\t ./up_count_lim <path_to_lh_proj> <proj_name> <cl>" << std::endl;
        return 1;
    }

    double cl;
    std::istringstream(argv[3]) >> cl;


    if(cl < 0 || cl > 1){
        std::cout << "Error: cl should be between 0 and 1" << std::endl;
        return 2;
    }

    // load up the root histogram
    TFile f(argv[1]);
    TH1D *rtHist = (TH1D*)f.Get(argv[2]);

    if(rtHist == NULL){
        std::cout << "Couldn't read TH1D " << argv[2] << " from " << argv[1]
                  << "\n\n Contents are: \n" << std::endl;
        f.ls();
        return 3;
    }
        

    Histogram lhProj  = DistTools::ToHist(*rtHist);

    double lim = BayesIntervalCalc::UpperBound(lhProj, cl);

    std::cout << "Upper limit is : " << lim << " counts @ " << 100 * cl 
              << "% " << std::endl;

    return 0;
}
