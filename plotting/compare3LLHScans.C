#include <iostream>
#include <filesystem>
#include <string>

/* ///////////////////////////////////////////////////////////////////
///
/// Script for comparing 3 LLH scans.
///
/// The user inputs the root files made by fixedOscLLHScan for
/// along with labels to be printed in the legend for
/// each.
///
/// Histograms in the input files are looped through and drawn.
///
/// The plots are drawn and saved in a pdf and root
/// file.
///
/////////////////////////////////////////////////////////////////// */

void LoopHistos(TDirectory *d1, TDirectory *d2, TDirectory *d3, std::string label1, std::string label2, std::string label3, std::string outfilename)
{

  // Get the directory
  TIter next(d1->GetListOfKeys());
  TKey *key;
  std::string name;

  // Loop through all entries
  while ((key = (TKey *)next()))
  {

    // Get object name and type
    std::string name = std::string(key->GetName());
    std::string classname = std::string(key->GetClassName());

    std::cout << name << std::endl;

    // If we're a histogram, let's plot it
    if (classname == "TH1F" || classname == "TH1D")
    {

      TCanvas *c = new TCanvas("c", "c", 1500, 1080);
      // Aesthetics
      gStyle->SetOptStat(0);
      c->SetTopMargin(0.1);
      c->SetBottomMargin(0.18);
      c->SetLeftMargin(0.18);
      c->SetRightMargin(0.09);
      c->SetGrid();
      c->SetFrameLineWidth(2);

      TH1D *plot1 = (TH1D *)d1->Get(name.c_str())->Clone();
      plot1->GetYaxis()->SetTitleOffset(1.5);
      c->cd();
      plot1->SetLineColor(kBlack);
      plot1->SetLineWidth(2);
      plot1->GetXaxis()->SetTitleSize(0.055);
      plot1->GetYaxis()->SetTitleSize(0.055);
      plot1->GetXaxis()->SetLabelSize(0.045);
      plot1->GetYaxis()->SetLabelSize(0.045);
      plot1->Draw();
      gPad->Update();

      TObject *obj2 = d2->Get(name.c_str());
      TH1D *plot2 = nullptr;
      if (obj2)
      {
        TH1D *hist2 = dynamic_cast<TH1D *>(obj2);
        if (hist2)
        {
          plot2->GetYaxis()->SetTitleOffset(1.5);
          c->cd();
          plot2->SetLineWidth(2);
          plot2->SetLineColor(kRed);
          plot2->SetLineStyle(2);
          plot2->GetXaxis()->SetTitleSize(0.055);
          plot2->GetYaxis()->SetTitleSize(0.055);
          plot2->GetXaxis()->SetLabelSize(0.045);
          plot2->GetYaxis()->SetLabelSize(0.045);
          plot2->Draw("same");
          gPad->Update();
        }
      }

      TObject *obj3 = d3->Get(name.c_str());
      TH1D *plot3 = nullptr;
      if (obj3)
      {
        TH1D *hist3 = dynamic_cast<TH1D *>(obj3);
        if (hist3)
        {
          plot3->GetYaxis()->SetTitleOffset(1.5);
          c->cd();
          plot3->SetLineWidth(2);
          plot3->SetLineStyle(3);
          plot3->SetLineColor(kBlue);
          plot3->GetXaxis()->SetTitleSize(0.055);
          plot3->GetYaxis()->SetTitleSize(0.055);
          plot3->GetXaxis()->SetLabelSize(0.045);
          plot3->GetYaxis()->SetLabelSize(0.045);
          plot3->Draw("same");
          gPad->Update();
        }
      }

      TLegend *t1 = new TLegend(0.4, 0.7, 0.6, 0.85);
      t1->SetLineWidth(2);
      t1->AddEntry(plot1, label1.c_str(), "l");
      t1->AddEntry(plot2, label2.c_str(), "l");
      t1->AddEntry(plot3, label3.c_str(), "l");
      t1->Draw("same");

      c->Print(outfilename.c_str());
      delete plot1;
      delete plot2;
      delete plot3;
      delete c;
    }
  }
}

void compare3LLHScans(std::string filename1, std::string filename2, std::string filename3, std::string label1, std::string label2, std::string label3)
{

  // Open file
  std::filesystem::path filepath1(filename1);
  TFile *f1 = new TFile(filename1.c_str(), "OPEN");
  TFile *f2 = new TFile(filename2.c_str(), "OPEN");
  TFile *f3 = new TFile(filename3.c_str(), "OPEN");
  std::string outputfilename = filepath1.replace_extension("comp.pdf").string();

  // Aesthetics
  gStyle->SetOptStat(0);

  // Print file name to first page
  TCanvas *c1 = new TCanvas("c1", "c1", 1500, 1080);
  c1->Print((outputfilename + "[").c_str());
  c1->cd();

  // Now loop over all the histos and print each
  LoopHistos(f1, f2, f3, label1, label2, label3, outputfilename);

  std::cout << "out the loop" << std::endl;

  c1->Print((outputfilename + "]").c_str());
}
