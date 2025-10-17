#include <TFile.h>
#include <TTree.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TStyle.h>

/* ///////////////////////////////////////////////////////////////////
///
/// Script for comparing 2D contours of two fits.
///
/// The user inputs the root files made by plotFixedOscLLH for
/// both fits, along with labels to be printed in the legend for
/// each, and an output filename.
///
/// The plot is drawn and the canvas is saved as a pdf and a root
/// file, along with each of the histograms.
///
/////////////////////////////////////////////////////////////////// */

void compareContours(std::string filename1, std::string label1, std::string filename2, std::string label2, std::string outfilename)
{

    // Open the ROOT file
    TFile *file1 = TFile::Open(filename1.c_str(), "READ");
    if (!file1 || file1->IsZombie())
    {
        std::cerr << "Error: Could not open file " << filename1 << std::endl;
        return;
    }

    TFile *file2 = TFile::Open(filename2.c_str(), "READ");
    if (!file2 || file2->IsZombie())
    {
        std::cerr << "Error: Could not open file " << filename2 << std::endl;
        return;
    }

    // Get the histos
    TH2D *hLLH1 = (TH2D *)file1->Get("hLLH");
    hLLH1->SetName("hLLH1");
    hLLH1->SetTitle("hLLH1");
    if (!hLLH1)
    {
        std::cerr << "No histogram hLLH found in file: " << filename1 << std::endl;
        file1->Close();
        delete file1;
        throw;
    }

    TH2D *hLLH2 = (TH2D *)file2->Get("hLLH");
    hLLH2->SetName("hLLH2");
    hLLH2->SetTitle("hLLH2");
    if (!hLLH2)
    {
        std::cerr << "No histogram hLLH found in file: " << filename2 << std::endl;
        file2->Close();
        delete file2;
        throw;
    }

    double contours[1];
    contours[0] = 2.295748928898636;

    hLLH1->SetLineColor(kBlue+2);
    hLLH1->SetLineWidth(2);
    hLLH1->SetLineStyle(1);
    hLLH1->SetContour(1, contours);
  
    hLLH2->SetLineColor(kRed+2);
    hLLH2->SetLineWidth(2);
    hLLH2->SetLineStyle(2);
    hLLH2->SetContour(1, contours);

    // Draw the histograms
    TCanvas *c1 = new TCanvas("c1", "Contours", 800, 600);
    c1->SetBottomMargin(0.13);
    gPad->SetFrameLineWidth(2);
    gStyle->SetOptStat(0);
    gPad->SetGrid(1);

    hLLH1->SetTitle("");
    hLLH1->Draw("cont2");
    hLLH2->Draw("cont3 SAME");

    TLegend *t1 = new TLegend(0.6, 0.75, 0.89, 0.87);
    t1->AddEntry(hLLH1, label1.c_str(), "l");
    t1->AddEntry(hLLH2, label2.c_str(), "l");
    t1->SetLineWidth(2);
    t1->SetTextFont(42);
    t1->Draw();

    // Save plot as image and rootfile
    c1->SaveAs((outfilename + ".pdf").c_str());
    c1->SaveAs((outfilename + ".root").c_str());
}
