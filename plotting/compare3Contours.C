#include <TFile.h>
#include <TTree.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TStyle.h>

/* ///////////////////////////////////////////////////////////////////
///
/// Script for comparing 2D contours of three fits.
///
/// The user inputs the root files made by plotFixedOscLLH for
/// each fit, along with labels to be printed in the legend for
/// each, and an output filename.
///
/// The plot is drawn and the canvas is saved as a pdf and a root
/// file, along with each of the histograms.
///
/////////////////////////////////////////////////////////////////// */

void compare3Contours(std::string filename1, std::string label1, std::string filename2, std::string label2,
                     std::string filename3, std::string label3, std::string outfilename)
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
    
    TFile *file3 = TFile::Open(filename3.c_str(), "READ");
    if (!file3 || file3->IsZombie())
    {
        std::cerr << "Error: Could not open file " << filename3 << std::endl;
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

    TH2D *hLLH3 = (TH2D *)file3->Get("hLLH");
    hLLH3->SetName("hLLH3");
    hLLH3->SetTitle("hLLH3");
    if (!hLLH3)
    {
        std::cerr << "No histogram hLLH found in file: " << filename3 << std::endl;
        file3->Close();
        delete file3;
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
    hLLH2->SetLineStyle(1);
    hLLH2->SetContour(1, contours);

    hLLH3->SetLineColor(kGreen+2);
    hLLH3->SetLineWidth(2);
    hLLH3->SetLineStyle(1);
    hLLH3->SetContour(1, contours);

    // Draw the histograms
    TCanvas *c1 = new TCanvas("c1", "Contours", 800, 600);
    c1->SetBottomMargin(0.13);
    gPad->SetFrameLineWidth(2);
    gStyle->SetOptStat(0);
    gPad->SetGrid(1);

    hLLH1->SetTitle("");
    hLLH1->Draw("cont2");
    hLLH2->Draw("cont3 SAME");
    hLLH3->Draw("cont3 SAME");

    TLegend *t1 = new TLegend(0.165, 0.17, 0.83, 0.34);
    t1->AddEntry(hLLH1, label1.c_str(), "l");
    t1->AddEntry(hLLH2, label2.c_str(), "l");
    t1->AddEntry(hLLH3, label3.c_str(), "l");
    t1->SetLineWidth(2);
    t1->SetTextFont(42);
    t1->SetTextSize(0.03);
    t1->Draw();

    // Save plot as image and rootfile
    c1->SaveAs((outfilename + ".pdf").c_str());
    c1->SaveAs((outfilename + ".root").c_str());
}
