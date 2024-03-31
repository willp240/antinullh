// compile me using $(source oxsx/bin/compile.sh to_h5.cpp)
// run with ./to_h5 <root-file> <new-h5-file>

#include <IO.h>
#include <DistTools.h>
#include <Histogram.h>
#include <TFile.h>
#include <TH2D.h>
#include <iostream>
#include <ContainerTools.hpp>

void convert(const std::string& infile_, const std::string& name_, const std::string& outfile_){
    TFile f(infile_.c_str());
    TH2D* h = (TH2D*)f.Get(name_.c_str());
    Histogram hist = DistTools::ToHist(*h);
    IO::SaveHistogram(hist, outfile_);
}
int main(int argc, char* argv[]){
    convert(argv[1], "", argv[2]);
    return 0;
}
