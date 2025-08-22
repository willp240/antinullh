#include <TFile.h>
#include <TTree.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TStyle.h>

/* ///////////////////////////////////////////////////////////////////
///
/// Script for comparing post fit parameter values and prefit
/// constraints relative to nominal values for three different sets of
/// fixed oscillation fits.
///
/// The user inputs the root files made by plotFixedOscParams for
/// all fits, along with labels to be printed in the legend for
/// each, and an output filename.
///
/// The nominal and constraint histograms are got from the first fit
/// params file, and postfit histograms are taken from each fit
/// param files.
///
/// The plot is drawn and the canvas is saved as a pdf and a root
/// file, along with each of the histograms.
///
/////////////////////////////////////////////////////////////////// */

void compare3FixedOscParams(std::string filename1, std::string label1, std::string filename2, std::string label2, std::string filename3, std::string label3, std::string outfilename)
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
    TH1D *hNom1 = (TH1D *)file1->Get("nominal");
    hNom1->SetName("nominal1");
    hNom1->SetTitle("nominal1");
    TH1D *hConstr1 = (TH1D *)file1->Get("constraints");
    hConstr1->SetName("constraints1");
    hConstr1->SetTitle("constraints1");
    TH1D *hPostFit1 = (TH1D *)file1->Get("postfit");
    hPostFit1->SetName("postfit1");
    hPostFit1->SetTitle("postfit1");

    TH1D *hPostFit2 = (TH1D *)file2->Get("postfit");
    hPostFit2->SetName("postfit2");
    hPostFit2->SetTitle("postfit2");

    TH1D *hPostFit3 = (TH1D *)file3->Get("postfit");
    hPostFit3->SetName("postfit3");
    hPostFit3->SetTitle("postfit3");

    // Draw the histograms
    TCanvas *c1 = new TCanvas("c1", "Params", 800, 600);
    c1->SetBottomMargin(0.18);
    gPad->SetFrameLineWidth(2);
    gStyle->SetOptStat(0);
    gPad->SetGrid(1);

    hNom1->SetLineColor(kGreen + 2);
    hNom1->SetLineWidth(12);
    hNom1->SetMarkerStyle(108);
    hNom1->SetMarkerColor(kGreen + 2);

    hConstr1->SetLineColor(kRed + 2);
    hConstr1->SetLineWidth(12);

    hPostFit1->SetLineColor(kBlue + 2);
    hPostFit1->SetLineWidth(8);
    hPostFit1->SetLineStyle(1);

    hPostFit2->SetLineColor(kGreen + 3);
    hPostFit2->SetLineWidth(4);
    hPostFit2->SetLineStyle(1);

    hPostFit3->SetLineColor(kBlack);
    hPostFit3->SetLineWidth(1);
    hPostFit3->SetLineStyle(1);

    hConstr1->GetYaxis()->SetRangeUser(0, 2);
    hConstr1->GetYaxis()->SetTitle("Relative to Nominal");
    hConstr1->GetYaxis()->SetTitleOffset(1.2);
    hConstr1->GetXaxis()->SetLabelOffset(0.007);
    hConstr1->GetXaxis()->SetTitle("Fit Parameters");
    hConstr1->GetXaxis()->SetTitleOffset(2.0);
    hConstr1->SetTitle("");

    hConstr1->GetXaxis()->SetTitleFont(42);
    hConstr1->GetYaxis()->SetTitleFont(42);
    hConstr1->GetXaxis()->SetLabelFont(42);
    hConstr1->GetYaxis()->SetLabelFont(42);
    hConstr1->SetTitleFont(42);

    hConstr1->Draw("E1");
    hPostFit1->Draw("E1same");
    hPostFit2->Draw("E1same");
    hPostFit3->Draw("E1same");

    TLegend *t1 = new TLegend(0.55, 0.63, 0.88, 0.88);
    t1->AddEntry(hConstr1, "Prefit", "lep");
    t1->AddEntry(hPostFit1, ("Postfit " + label1).c_str(), "lep");
    t1->AddEntry(hPostFit2, ("Postfit " + label2).c_str(), "lep");
    t1->AddEntry(hPostFit3, ("Postfit " + label3).c_str(), "lep");
    t1->SetLineWidth(2);
    t1->SetTextFont(42);
    t1->Draw();

    // Save plot as image and rootfile
    c1->SaveAs((outfilename + ".pdf").c_str());
    c1->SaveAs((outfilename + ".root").c_str());
}
