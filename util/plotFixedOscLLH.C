#include <TFile.h>
#include <TTree.h>
#include <TH2D.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <iostream>

std::pair<TH1D *, TH1D *> GetMinProfiles(TH2D *hist2D)
{

    std::pair<TH1D *, TH1D *> Profile;

    int nBinsX = hist2D->GetNbinsX();
    int nBinsY = hist2D->GetNbinsY();

    // Create a TH1D histogram to store the minimum Y value for each X bin
    TH1D *minProfileX = new TH1D("minProfileX", "minProfileX", nBinsX, hist2D->GetXaxis()->GetXmin(), hist2D->GetXaxis()->GetXmax());
    TH1D *minProfileY = new TH1D("minProfileY", "minProfileY", nBinsX, hist2D->GetYaxis()->GetXmin(), hist2D->GetYaxis()->GetXmax());
    std::cout << hist2D->GetXaxis()->GetXmin() << " " << hist2D->GetXaxis()->GetXmax() << std::endl;
    std::cout << hist2D->GetYaxis()->GetXmin() << " " << hist2D->GetYaxis()->GetXmax() << std::endl;

    for (int i = 1; i <= nBinsX; i++)
    {
        double minVal = std::numeric_limits<double>::max(); // Set to large value
        bool hasValue = false;                              // To check if any non-zero bin exists

        // Loop over all Y bins at this X
        for (int j = 1; j <= nBinsY; j++)
        {
            double binContent = hist2D->GetBinContent(i, j);
            if (binContent != 0 && binContent < minVal)
            {
                minVal = binContent;
                hasValue = true;
            }
        }

        // Set the min value in the new histogram (use 0 if no value found)
        if (hasValue)
        {
            minProfileX->SetBinContent(i, minVal);
        }
        else
        {
            minProfileX->SetBinContent(i, 0);
        }
    }
    Profile.first = minProfileX;
    minProfileX->SetLineWidth(2);

    // Loop over each Y bin
    for (int i = 1; i <= nBinsY; i++)
    {
        double minVal = std::numeric_limits<double>::max(); // Set to large value
        bool hasValue = false;                              // To check if any non-zero bin exists

        // Loop over all Y bins at this X
        for (int j = 1; j <= nBinsX; j++)
        {
            double binContent = hist2D->GetBinContent(j, i);
            if (binContent != 0 && binContent < minVal)
            {
                minVal = binContent;
                hasValue = true;
            }
        }

        // Set the min value in the new histogram (use 0 if no value found)
        if (hasValue)
            minProfileY->SetBinContent(i, minVal);
        else
            minProfileY->SetBinContent(i, 0);
    }
    minProfileY->SetLineWidth(2);
    Profile.second = minProfileY;

    return Profile;
}

std::map<std::string, double> extractXForMinYAndYDiff_1sigma(TH1D *hist1D)
{
    std::cout << "inside " << std::endl;

    double minY = std::numeric_limits<double>::max(); // hist1D->GetMinimum();
    int min_index = -1;

    double minX = -1; // hist1D->GetBinCenter(hist1D->GetMinimumBin());
    std::map<std::string, double> x_map;
    // Get the number of bins in the histogram
    int nBins = hist1D->GetNbinsX();

    for (int i = 1; i < nBins; ++i)
    { // Bins start at 1 in ROOT
        double yValue = hist1D->GetBinContent(i);
        //std::cout << i << " " << yValue << " " << minY << " " << min_index << std::endl;
        // std::cout<<hist1D->GetBinCenter(i)<<" "<<yValue<<std::endl;
        if (yValue < minY)
        {
            minY = yValue;
            min_index = i; // Store the bin index of the minimum
        }
    }
    minX = hist1D->GetBinCenter(min_index);
    std::cout << "minY " << minY << std::endl;
    std::cout << "min_index " << min_index << std::endl;
    std::cout << "minX " << minX << std::endl;
    x_map["nom"] = minX;

    // Find Right boundary: looping from min_index to the end
    for (int i = min_index + 1; i < nBins; i++)
    {
        double yCurrent = hist1D->GetBinContent(i);
        double yNext = hist1D->GetBinContent(i + 1);
        double xCurrent = hist1D->GetBinCenter(i);
        double xNext = hist1D->GetBinCenter(i + 1);
        if (std::abs(yCurrent - minY) <= 1 && std::abs(yNext - minY) >= 1)
        {
            x_map["right"] = xNext;
            break;
        }
    }
    // Find left boundary: looping from 0 to min_index
    for (int i = 1; i < min_index; i++)
    {
        double yCurrent = hist1D->GetBinContent(i);
        double yNext = hist1D->GetBinContent(i + 1);
        double xCurrent = hist1D->GetBinCenter(i);
        double xNext = hist1D->GetBinCenter(i + 1);

        if (std::abs(yCurrent - minY) >= 1 && std::abs(yNext - minY) <= 1)
        {
            x_map["left"] = xCurrent;
        }
    }
    std::cout << x_map["right"] << " " << x_map["left"] << std::endl;

    // Return the minimum X value and the pairs of X values where Y - Ymin = 1
    return x_map;
}

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

    hLLH->GetXaxis()->SetTitleFont(42);
    hLLH->GetYaxis()->SetTitleFont(42);
    hLLH->GetZaxis()->SetTitleFont(42);
    hLLH->GetXaxis()->SetLabelFont(42);
    hLLH->GetYaxis()->SetLabelFont(42);
    hLLH->GetZaxis()->SetLabelFont(42);
    hLLH->SetTitleFont(42);

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

    std::pair<TH1D *, TH1D *> Profile = GetMinProfiles(hLLH);

    // Profile LLH of theta
    TCanvas *c2 = new TCanvas("c2", "theta LLH", 800, 600);
    c2->SetRightMargin(0.15);
    gPad->SetFrameLineWidth(2);
    gStyle->SetOptStat(0);
    gPad->SetGrid();

    Profile.first->Draw();
    Profile.first->SetLineColor(kBlack);
    Profile.first->GetXaxis()->SetTitle("#theta_{12}");
    Profile.first->GetYaxis()->SetTitle("2#Deltaln(L)");
    Profile.first->GetYaxis()->SetTitleOffset(1.2);
    Profile.first->SetTitle("");

    std::map<std::string, double> DelTheta_XVal = extractXForMinYAndYDiff_1sigma(Profile.first);
    std::cout << DelTheta_XVal["nom"] << " " << DelTheta_XVal["left"] << " " << DelTheta_XVal["right"] << std::endl;
    double yMin1 = 0;
    double yMax1 = Profile.first->GetMaximum();
    double yMaxCanvas1 = yMax1 * 1.2;
    Profile.first->SetMaximum(yMaxCanvas1);
    TLine *leftLine1 = new TLine(DelTheta_XVal["left"], yMin1, DelTheta_XVal["left"], yMaxCanvas1);
    TLine *rightLine1 = new TLine(DelTheta_XVal["right"], yMin1, DelTheta_XVal["right"], yMaxCanvas1);
    leftLine1->SetLineColor(kBlue);
    leftLine1->SetLineStyle(2);
    leftLine1->SetLineWidth(2);
    rightLine1->SetLineColor(kBlue);
    rightLine1->SetLineStyle(2);
    rightLine1->SetLineWidth(2);

    leftLine1->Draw("same");
    rightLine1->Draw("same");
    TLatex latex1;
    latex1.SetTextSize(0.035); // Adjust text size
    latex1.SetTextFont(42);
    latex1.SetTextColor(kBlack);
    double left_sigma1 = DelTheta_XVal["nom"] - DelTheta_XVal["left"];
    double right_sigma1 = DelTheta_XVal["right"] - DelTheta_XVal["nom"];
    latex1.DrawLatex(DelTheta_XVal["right"] - 30, yMaxCanvas1 * 0.9, Form("#theta_{12} = %.2f#circ ^{+%.2f}_{-%.2f}", DelTheta_XVal["nom"], right_sigma1, left_sigma1));
    
    pathObj.replace_filename("theta12LLHDiff.pdf");
    c2->SaveAs(pathObj.string().c_str());
    pathObj.replace_filename("theta12LLHDiff.root");
    c2->SaveAs(pathObj.string().c_str());

    // Profile LLH on deltaM
    TCanvas *c3 = new TCanvas("c3", "DeltaM LLH", 800, 600);
    c3->SetRightMargin(0.15);
    gPad->SetFrameLineWidth(2);
    gStyle->SetOptStat(0);
    gPad->SetGrid();

    Profile.second->Draw();
    Profile.second->SetLineColor(kBlack);
    Profile.second->GetXaxis()->SetTitle("#Delta m^{2}_{21}, 10^{-5}MeV");
    Profile.second->GetYaxis()->SetTitle("2#Deltaln(L)");
    Profile.second->GetYaxis()->SetTitleOffset(1.2);
    Profile.second->SetTitle("");

    std::map<std::string, double> DelM_XVal = extractXForMinYAndYDiff_1sigma(Profile.second);
    double yMin = 0;
    double yMax = Profile.second->GetMaximum();
    double yMaxCanvas = yMax * 1.2;
    Profile.second->SetMaximum(yMaxCanvas);
    TLine *leftLine = new TLine(DelM_XVal["left"], yMin, DelM_XVal["left"], yMaxCanvas);
    TLine *rightLine = new TLine(DelM_XVal["right"], yMin, DelM_XVal["right"], yMaxCanvas);
    leftLine->SetLineColor(kBlue);
    leftLine->SetLineStyle(2);
    leftLine->SetLineWidth(2);
    rightLine->SetLineColor(kBlue);
    rightLine->SetLineStyle(2);
    rightLine->SetLineWidth(2);

    leftLine->Draw("same");
    rightLine->Draw("same");
    TLatex latex;
    latex.SetTextSize(0.035); // Adjust text size
    latex.SetTextFont(42);
    latex.SetTextColor(kBlack);
    double left_sigma = DelM_XVal["nom"] - DelM_XVal["left"];
    double right_sigma = DelM_XVal["right"] - DelM_XVal["nom"];
    latex.DrawLatex(DelM_XVal["left"] - 0.000033, yMaxCanvas * 0.9, Form("#Deltam^{2}_{21} = %.3f^{+%.3f}_{-%.3f}#times10^{-5} eV^{2}", DelM_XVal["nom"] * 1E5, right_sigma * 1E5, left_sigma * 1E5));
    pathObj.replace_filename("delM12LLHDiff.pdf");
    c3->SaveAs(pathObj.string().c_str());
    pathObj.replace_filename("delM12LLHDiff.root");
    c3->SaveAs(pathObj.string().c_str());

    file->Close();
}
