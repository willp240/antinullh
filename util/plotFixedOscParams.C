#include <TFile.h>
#include <TTree.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TStyle.h>

void plotFixedOscParams(const char *filename = "fit_results.root")
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
    file->GetObject("FitResults", tree);
    if (!tree)
    {
        std::cerr << "Error: Could not find TTree 'fit_tree' in file " << filename << std::endl;
        file->Close();
        return;
    }

    // Define variables to hold the tree branches
    double theta, deltam, llh;
    tree->SetBranchAddress("theta", &theta);
    tree->SetBranchAddress("deltam", &deltam);
    tree->SetBranchAddress("LLH", &llh);

    int nEntries = tree->GetEntries();

    // Create a TH2D histogram
    int nBinsX = sqrt(nEntries); // Adjust binning as needed
    int nBinsY = sqrt(nEntries);
    double minDeltam = tree->GetMinimum("deltam");
    double maxDeltam = tree->GetMaximum("deltam");
    double minTheta = tree->GetMinimum("theta");
    double maxTheta = tree->GetMaximum("theta");
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
