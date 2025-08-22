// ROOT Headers
#include <TFile.h>
#include <TTree.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TStyle.h>

// c++ Headers
#include <sys/stat.h>

/* ///////////////////////////////////////////////////////////////////
///
/// Script for plotting the unscaled PDFs inside a directory,
/// projected onto the Y axis (alpha-n classifier). The user
/// inputs the directory name, and it gets the root files in that
/// directory and plots the 1D projections onto one plot, which is
/// saved as a PDF
///
/////////////////////////////////////////////////////////////////// */

void projectPDFs(std::string dirname)
{
    gStyle->SetOptStat(0);

    std::vector<std::filesystem::path> pdfFiles;
    std::string outputfilename = dirname + "/pdfplots.pdf";

    // Print file name to first page
    TCanvas *c1 = new TCanvas("c1", "c1", 1500, 1080);
    c1->SetTopMargin(0.1);
    c1->SetBottomMargin(0.18);
    c1->SetLeftMargin(0.15);
    c1->SetRightMargin(0.07);
    c1->SetGrid();
    c1->SetFrameLineWidth(2);
    c1->cd();

    TFile *fAlphaC = TFile::Open((dirname + "/alphan_CScatter.root").c_str(), "READ");
    TH2D *hAlphaC = (TH2D *)fAlphaC->Get("alphan_CScatter");
    hAlphaC->SetLineColor(kBlue);
    hAlphaC->ProjectionY()->SetLineWidth(2);
    hAlphaC->ProjectionY()->SetMaximum(0.1);
    hAlphaC->ProjectionY()->GetYaxis()->SetTitle("Normalised Probability");
    hAlphaC->ProjectionY()->GetYaxis()->SetTitleFont(42);
    hAlphaC->ProjectionY()->GetYaxis()->SetLabelFont(42);
    hAlphaC->ProjectionY()->GetYaxis()->SetTitleOffset(1.2);
    hAlphaC->ProjectionY()->Draw("hist");
    TFile *fAlphaO = TFile::Open((dirname + "/alphan_OExcited.root").c_str(), "READ");
    TH2D *hAlphaO = (TH2D *)fAlphaO->Get("alphan_OExcited");
    hAlphaO->SetLineColor(kBlue + 2);
    hAlphaO->ProjectionY()->SetLineWidth(2);
    hAlphaO->ProjectionY()->Draw("histsame");
    TFile *fAlphaP = TFile::Open((dirname + "/alphan_PRecoil.root").c_str(), "READ");
    TH2D *hAlphaP = (TH2D *)fAlphaP->Get("alphan_PRecoil");
    hAlphaP->SetLineColor(kAzure + 10);
    hAlphaP->ProjectionY()->SetLineWidth(2);
    hAlphaP->ProjectionY()->Draw("histsame");
    TFile *fGeoU = TFile::Open((dirname + "/geonu_U.root").c_str(), "READ");
    TH2D *hGeoU = (TH2D *)fGeoU->Get("geonu_U");
    hGeoU->SetLineColor(kGreen);
    hGeoU->ProjectionY()->SetLineWidth(2);
    hGeoU->ProjectionY()->Draw("histsame");
    TFile *fGeoTh = TFile::Open((dirname + "/geonu_Th.root").c_str(), "READ");
    TH2D *hGeoTh = (TH2D *)fGeoTh->Get("geonu_Th");
    hGeoTh->SetLineColor(kGreen + 2);
    hGeoTh->ProjectionY()->SetLineWidth(2);
    hGeoTh->ProjectionY()->Draw("histsame");
    TFile *fReac = TFile::Open((dirname + "/reactor_nubar.root").c_str(), "READ");
    TH2D *hReac = (TH2D *)fReac->Get("reactor_nubar");
    hReac->SetLineColor(kRed);
    hReac->ProjectionY()->SetLineWidth(2);
    hReac->ProjectionY()->Draw("histsame");

    TLegend *t1 = new TLegend(0.7, 0.55, 0.88, 0.85);
    t1->SetLineWidth(2);
    t1->AddEntry(hReac, "Reac-#bar{#nu}", "l");
    t1->AddEntry(hGeoU, "Geo-U", "l");
    t1->AddEntry(hGeoTh, "Geo-Th", "l");
    t1->AddEntry(hAlphaP, "#alpha-n P. Recoil", "l");
    t1->AddEntry(hAlphaC, "#alpha-n O Excitation", "l");
    t1->AddEntry(hAlphaO, "#alpha-n C Scattering", "l");
    t1->SetTextFont(42);
    t1->Draw();

    c1->SaveAs((dirname + "/projections.pdf").c_str());
    c1->SaveAs((dirname + "/projections.root").c_str());
}
