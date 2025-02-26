#include <TFile.h>
#include <TTree.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TStyle.h>

/* ///////////////////////////////////////////////////////////////////
///
/// Script for plotting bestfit LLH from fixed oscillation fits for
/// grid scan of oscillation parameters.
/// 
/// The user inputs the root file made by makeFixedOscTree, and it 
/// reads the oscillation parameters and LLH for each entry. These
/// are drawn in a TH2D, along with contours.asm
///
/// The plot is drawn and the canvas is saved as a root and pdf file
///
/////////////////////////////////////////////////////////////////// */

void plotFixedOscLLH(const char *filename = "fit_results.root")
{

    // Open the ROOT file
    TFile *file = TFile::Open(filename, "READ");
    if (!file || file->IsZombie())
    {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        return;
    }

    // Get the TTree
    TTree *tree = nullptr;
    file->GetObject("fitResults", tree);
    if (!tree)
    {
        std::cerr << "Error: Could not find TTree 'fit_tree' in file " << filename << std::endl;
        file->Close();
        return;
    }

    // Define variables to hold the tree branches
    double theta, deltam, llh;
    tree->SetBranchAddress("theta12", &theta);
    tree->SetBranchAddress("deltam21", &deltam);
    tree->SetBranchAddress("LLH", &llh);

    int nEntries = tree->GetEntries();

    // Create a TH2D histogram
    int nBinsX = sqrt(nEntries);
    int nBinsY = sqrt(nEntries);
    double minDeltam = tree->GetMinimum("deltam21");
    double maxDeltam = tree->GetMaximum("deltam21");
    double minTheta = tree->GetMinimum("theta12");
    double maxTheta = tree->GetMaximum("theta12");
    double minLLH = tree->GetMinimum("LLH");

    TH2D *hLLH = new TH2D("hLLH", "#Delta LLH;#Delta m^2, MeV;#theta",
                          nBinsX, minTheta, maxTheta,
                          nBinsY, minDeltam, maxDeltam);

    // Loop over the tree entries and fill the histogram
    for (Long64_t i = 0; i < nEntries; i++)
    {
        tree->GetEntry(i);
        hLLH->Fill(theta, deltam, 2 * (llh - minLLH));
    }

    // Draw the histogram
    TCanvas *c1 = new TCanvas("c1", "LLH", 800, 600);
    c1->SetRightMargin(0.15);
    gPad->SetFrameLineWidth(2);
    gStyle->SetOptStat(0);

    hLLH->GetXaxis()->SetTitle("#theta_{12}");
    hLLH->GetYaxis()->SetTitle("#Delta m^{2}_{21}, MeV");
    hLLH->GetYaxis()->SetTitleOffset(1.2);
    hLLH->GetZaxis()->SetTitle("2#Deltaln(L)");
    hLLH->SetTitle("");

    double contours[1];
    contours[0] = 2.295748928898636;

    hLLH->DrawCopy("colz");
    hLLH->SetContour(1, contours);
    hLLH->Draw("cont3 same");
    hLLH->SetLineColor(kRed);

    // Save plot as image and rootfile
    std::filesystem::path pathObj(filename);
    pathObj.replace_filename("LLH2D.pdf");
    c1->SaveAs(pathObj.string().c_str());
    pathObj.replace_filename("LLH2D.root");
    c1->SaveAs(pathObj.string().c_str());

    file->Close();
}
