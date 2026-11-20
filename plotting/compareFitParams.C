#include "FitParamDrawing.hh"

/* ///////////////////////////////////////////////////////////////////
///
/// Compare parameter values from two fixed-oscillation, floated Minuit or MCMC fits
/// Read params.root files produced by plotFitParams or plotFixedOscParams
/// Use nominal values and constraints from the first fit
/// Save the comparison canvas and histograms as PDF and ROOT files
///
/////////////////////////////////////////////////////////////////// */

void compareFitParams(std::string filename1, std::string label1, std::string filename2, std::string label2, std::string outfilename)
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

    // Draw the histograms
    TCanvas *c1 = new TCanvas("c1", "Params", 800, 600);
    c1->SetBottomMargin(0.18);
    c1->SetTopMargin(0.22);
    gPad->SetFrameLineWidth(2);
    gStyle->SetOptStat(0);
    gPad->SetGrid(1);

    hNom1->SetLineColor(kGreen);
    hNom1->SetLineWidth(2);
    hConstr1->SetLineColor(kBlue);
    hConstr1->SetLineWidth(2);
    hPostFit1->SetLineColor(kRed);
    hPostFit1->SetLineWidth(2);
    hPostFit2->SetLineColor(kBlack);
    hPostFit2->SetLineWidth(2);
    hPostFit2->SetLineStyle(2);

    TGraphAsymmErrors *intervals1 = antinupostfit::ReadFitIntervals(*file1, hPostFit1, hConstr1, "postfit1");
    intervals1->SetMarkerStyle(20);
    TGraphAsymmErrors *intervals2 = antinupostfit::ReadFitIntervals(*file2, hPostFit2, hConstr1, "postfit2");
    intervals2->SetMarkerStyle(21);
    TGraphAsymmErrors *regions1 = antinupostfit::ReadFitRegions(*file1, hPostFit1, "regions1");
    TGraphAsymmErrors *regions2 = antinupostfit::ReadFitRegions(*file2, hPostFit2, "regions2");
    antinupostfit::SetFitComparisonRange(*hConstr1, {intervals1, intervals2, regions1, regions2});
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
    if (regions1->GetN())
        regions1->Draw("Z SAME");
    if (regions2->GetN())
        regions2->Draw("Z SAME");
    intervals1->Draw("PZ SAME");
    intervals2->Draw("PZ SAME");

    TLegend *t1 = new TLegend(0.12, 0.80, 0.88, 0.97);
    t1->AddEntry(hConstr1, "Prefit", "l");
    t1->AddEntry(intervals1, ("Postfit " + label1 + antinupostfit::FitIntervalLabel(*file1)).c_str(), "l");
    t1->AddEntry(intervals2, ("Postfit " + label2 + antinupostfit::FitIntervalLabel(*file2)).c_str(), "l");
    t1->SetBorderSize(0);
    t1->SetFillStyle(0);
    t1->SetTextSize(0.025);
    t1->SetLineWidth(2);
    t1->SetTextFont(42);
    t1->Draw();

    // Save plot as image and rootfile
    c1->SaveAs((outfilename + ".pdf").c_str());
    c1->SaveAs((outfilename + ".root").c_str());
}
