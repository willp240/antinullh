#include <iostream>
#include <filesystem>
#include <string>

void LoopHistos(TDirectory *dir, std::string outfilename) {

  // Get the directory
  TIter next(dir->GetListOfKeys());
  TKey *key;
  std::string name;

  // Loop through all entries
  while ((key = (TKey*)next())) {
    
    // Get object name and type
    std::string name = std::string(key->GetName());
    std::string classname = std::string(key->GetClassName());

    std::cout << name << std::endl;
    
    // If we're in a directory, recurse
    if (classname == "TDirectoryFile") {
      dir->cd(name.c_str());
      TDirectory *SubDir = gDirectory;
      LoopHistos(SubDir, outfilename);
    } 

    // If we're a histogram, let's plot it
    if (classname == "TH1F" || classname == "TH1D") {
      
      TCanvas *c = new TCanvas("c", "c", 1500, 1080);
      // Aesthetics
      gStyle->SetOptStat(0);
      c->SetTopMargin(0.1);
      c->SetBottomMargin(0.08);
      c->SetLeftMargin(0.15);
      c->SetRightMargin(0.06);
      c->SetGrid();
      c->SetFrameLineWidth(2);
      
      TH1D *plot1 = (TH1D*)dir->Get(name.c_str())->Clone();
      plot1->GetYaxis()->SetTitleOffset(2.0);
      c->cd();
      plot1->SetLineWidth(2);
      plot1->Draw();
      gPad->Update();
      c->Print(outfilename.c_str());
      delete plot1;
      delete c;
    }
  }
}

void PlotLLH(std::string filename) {

  // Open file
  TFile* File = new TFile(filename.c_str(), "OPEN");
  std::filesystem::path filepath(filename);
  std::string outputfilename = filepath.replace_extension(".pdf").string();

  // Aesthetics
  gStyle->SetOptStat(0);
  
  // Print file name to first page
  TCanvas *c1 = new TCanvas("c1", "c1", 1500, 1080);
  c1->Print((outputfilename+"[").c_str());
  c1->cd();

  // Now loop over all the histos and print each
  LoopHistos(File, outputfilename);

  std::cout << "out the loop" << std::endl;
  
  c1->Print((outputfilename+"]").c_str());
}
