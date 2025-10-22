// ROOT Headers
#include <TFile.h>
#include <TTree.h>
#include <TH2D.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TStyle.h>

// c++ Headers
#include <iostream>
#include <sys/stat.h>

/* ///////////////////////////////////////////////////////////////////
///
/// Script for plotting bestfit LLH from fixed oscillation fits for
/// grid scan of oscillation parameters.
///
/// The user inputs the root file made by makeFixedOscTree, and it
/// reads the oscillation parameters and LLH for each entry. These
/// are drawn in a TH2D, along with contours.
///
/// Profile LLH are calculated and plotted for each axis, and the
/// assymetric 1 sigma bounds are calculated and saved in the
/// outputted root file.
///
/// The plots are drawn and the canvases are saved as a root and
/// pdf files
///
/////////////////////////////////////////////////////////////////// */

// Function that takes in a TH2D, and returns a pair of TH1Ds. Each is a profile of the two axis of the TH2D
// It loops through the bins on one axis. For each bin, it loops through all the bins on the other axis,
// finding the minimum LLH for that row. This is then repeated for the other axis
std::pair<TH1D *, TH1D *> GetMinProfiles(TH2D *hist2D)
{

    std::pair<TH1D *, TH1D *> Profile;

    int nBinsX = hist2D->GetNbinsX();
    int nBinsY = hist2D->GetNbinsY();

    // Create a TH1D histogram to store the minimum Y value for each X bin
    TH1D *minProfileX = new TH1D("minProfileX", "minProfileX", nBinsX, hist2D->GetXaxis()->GetXmin(), hist2D->GetXaxis()->GetXmax());
    TH1D *minProfileY = new TH1D("minProfileY", "minProfileY", nBinsX, hist2D->GetYaxis()->GetXmin(), hist2D->GetYaxis()->GetXmax());

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

// Function that takes in a TH1D, and returns the best LLH and assymetric errors. It first finds the bin that has lowest profiled LLH,
// and counts out from there to find where LLH changes by more than 1.0 on each side
std::map<std::string, double> calc1Sigma(TH1D *hist1D)
{

    // Here Y refers to the LLH axis, X refers to the parameter axis
    double minY = std::numeric_limits<double>::max();
    int min_index = -1;

    double minX = -1;
    std::map<std::string, double> x_map;

    // Get the number of bins in the histogram
    int nBins = hist1D->GetNbinsX();

    for (int i = 1; i < nBins; ++i)
    {
        double yValue = hist1D->GetBinContent(i);
        if (yValue < minY)
        {
            minY = yValue;
            min_index = i; // Store the bin index of the minimum
        }
    }

    minX = hist1D->GetBinCenter(min_index);
    x_map["bestfit"] = minX;

    // Find right boundary: looping from min_index to the end
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

    // Return the minimum X value and the pairs of X values where Y - Ymin = 1
    return x_map;
}

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

    // Check the form of theta 12
    std::string theta12name = "";
    std::string theta12label = "";
    std::string theta12unit = "";
    std::string theta12labelunit = "";
    if (tree->GetBranch("theta12"))
    {
        theta12name = "theta12";
        theta12label = "#theta_{12}";
        theta12unit = "#circ";
        theta12labelunit = theta12label + ", " + theta12unit;
    }
    else if (tree->GetBranch("sintheta12"))
    {
        theta12name = "sintheta12";
        theta12label = "sin#theta_{12}";
        theta12labelunit = theta12label;
    }
    else if (tree->GetBranch("sinsqtheta12"))
    {
        theta12name = "sinsqtheta12";
        theta12label = "sin^{2}#theta_{12}";
        theta12labelunit = theta12label;
    }

    // Define variables to hold the tree branches
    double theta, deltam, llh;
    tree->SetBranchAddress(theta12name.c_str(), &theta);
    tree->SetBranchAddress("deltam21", &deltam);
    tree->SetBranchAddress("LLH", &llh);

    int nEntries = tree->GetEntries();

    // Create a TH2D histogram
    int nBinsX = sqrt(nEntries);
    int nBinsY = sqrt(nEntries);
    double minDeltam = tree->GetMinimum("deltam21");
    double maxDeltam = tree->GetMaximum("deltam21");
    double minTheta = tree->GetMinimum(theta12name.c_str());
    double maxTheta = tree->GetMaximum(theta12name.c_str());
    double minLLH = tree->GetMinimum("LLH");
    
    // Steps between centers (because both min and max are included as centers):
    double stepTheta = (maxTheta   - minTheta)   / (nBinsX - 1);
    double stepDm    = (maxDeltam  - minDeltam)  / (nBinsY - 1);

    double minThetaBin = minTheta-(stepTheta/2);
    double maxThetaBin = maxTheta+(stepTheta/2);
    double minDeltamBin = minDeltam-(stepDm/2);
    double maxDeltamBin = maxDeltam+(stepDm/2);

    // Make the histogram with those edges:
    TH2D* hLLH = new TH2D("hLLH", ("#Delta LLH;#Delta m^{2}, eV^{2};" + theta12labelunit).c_str(),
			  nBinsX, minThetaBin, maxThetaBin,
			  nBinsY, minDeltamBin, maxDeltamBin);
    
    // Loop over the tree entries and fill the histogram
    for (Long64_t i = 0; i < nEntries; i++)
    {
        tree->GetEntry(i);

	    int binx = hLLH->GetXaxis()->FindBin(theta);
	    int biny = hLLH->GetYaxis()->FindBin(deltam);

	    hLLH->SetBinContent(binx, biny, 2*(llh - minLLH));
    }
    
    // Draw the 2D histogram
    TCanvas *c1 = new TCanvas("c1", "LLH", 800, 600);
    c1->SetRightMargin(0.15);
    gPad->SetFrameLineWidth(2);
    gStyle->SetOptStat(0);

    hLLH->GetXaxis()->SetTitle(theta12labelunit.c_str());
    hLLH->GetYaxis()->SetTitle("#Delta m^{2}_{21}, eV^{2}");
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

    // Save plot as image and rootfile
    struct stat st = {0};
    std::filesystem::path pathObj(filename);
    pathObj.replace_filename("plots/");
    if (stat(pathObj.string().c_str(), &st) == -1)
        mkdir(pathObj.string().c_str(), 0700);
    pathObj.replace_filename("LLH.root");
    TFile *outfile = new TFile(pathObj.string().c_str(), "RECREATE");
    outfile->cd();
    hLLH->Write("hLLH");

    hLLH->DrawCopy("colz");
    hLLH->SetContour(1, contours);
    hLLH->Draw("cont3 same");
    hLLH->SetLineColor(kRed);
    pathObj.replace_filename("LLH2D.pdf");
    c1->SaveAs(pathObj.string().c_str());
    hLLH->Write("hLLH");
    c1->Write("c1");

    // Now we're going to make the profile LLH plots
    std::pair<TH1D *, TH1D *> Profile = GetMinProfiles(hLLH);

    // Profile LLH of theta
    TCanvas *c2 = new TCanvas("c2", "theta LLH", 800, 600);
    c2->SetRightMargin(0.15);
    gPad->SetFrameLineWidth(2);
    gStyle->SetOptStat(0);
    gPad->SetGrid();

    Profile.first->Draw();
    Profile.first->SetLineColor(kBlack);
    Profile.first->GetXaxis()->SetTitle(theta12labelunit.c_str());
    Profile.first->GetYaxis()->SetTitle("2#Deltaln(L)");
    Profile.first->GetYaxis()->SetTitleOffset(1.2);
    Profile.first->SetTitle("");
    Profile.first->GetXaxis()->SetTitleFont(42);
    Profile.first->GetYaxis()->SetTitleFont(42);
    Profile.first->GetZaxis()->SetTitleFont(42);
    Profile.first->GetXaxis()->SetLabelFont(42);
    Profile.first->GetYaxis()->SetLabelFont(42);
    Profile.first->GetZaxis()->SetLabelFont(42);
    Profile.first->SetTitleFont(42);

    std::map<std::string, double> thSigmas = calc1Sigma(Profile.first);
    double yMin = 0;
    double yMax = Profile.first->GetMaximum();
    double yMaxCanvas = yMax * 1.2;
    Profile.first->SetMaximum(yMaxCanvas);
    TLine *leftLine = new TLine(thSigmas["left"], yMin, thSigmas["left"], yMaxCanvas);
    TLine *rightLine = new TLine(thSigmas["right"], yMin, thSigmas["right"], yMaxCanvas);
    leftLine->SetLineColor(kBlue);
    leftLine->SetLineStyle(2);
    leftLine->SetLineWidth(2);
    rightLine->SetLineColor(kBlue);
    rightLine->SetLineStyle(2);
    rightLine->SetLineWidth(2);

    leftLine->Draw("same");
    rightLine->Draw("same");
    TLatex latex;
    latex.SetTextSize(0.035);
    latex.SetTextFont(42);
    latex.SetTextColor(kBlack);
    double left_sigma = thSigmas["bestfit"] - thSigmas["left"];
    double right_sigma = thSigmas["right"] - thSigmas["bestfit"];

    double labelpos = 0;
    if (theta12name == "theta12")
        labelpos = 35;
    else
        labelpos = 0.42;

    latex.DrawLatex(labelpos, yMaxCanvas * 0.9, Form("%s = (%.2f ^{+%.2f}_{-%.2f}) %s", theta12label.c_str(), thSigmas["bestfit"], right_sigma, left_sigma, theta12unit.c_str()));
    pathObj.replace_filename("theta12LLHDiff.pdf");
    c2->SaveAs(pathObj.string().c_str());
    outfile->cd();
    c2->Write("c2");
    Profile.first->Write("hTheta");

    // Profile LLH on deltaM
    TCanvas *c3 = new TCanvas("c3", "DeltaM LLH", 800, 600);
    c3->SetRightMargin(0.15);
    gPad->SetFrameLineWidth(2);
    gStyle->SetOptStat(0);
    gPad->SetGrid();

    Profile.second->Draw();
    Profile.second->SetLineColor(kBlack);
    Profile.second->GetXaxis()->SetTitle("#Delta m^{2}_{21}, 10^{-5}eV^{2}");
    Profile.second->GetYaxis()->SetTitle("2#Deltaln(L)");
    Profile.second->GetYaxis()->SetTitleOffset(1.2);
    Profile.second->SetTitle("");
    Profile.second->GetXaxis()->SetTitleFont(42);
    Profile.second->GetYaxis()->SetTitleFont(42);
    Profile.second->GetZaxis()->SetTitleFont(42);
    Profile.second->GetXaxis()->SetLabelFont(42);
    Profile.second->GetYaxis()->SetLabelFont(42);
    Profile.second->GetZaxis()->SetLabelFont(42);
    Profile.second->SetTitleFont(42);

    std::map<std::string, double> dmSigmas = calc1Sigma(Profile.second);
    yMin = 0;
    yMax = Profile.second->GetMaximum();
    yMaxCanvas = yMax * 1.2;
    Profile.second->SetMaximum(yMaxCanvas);
    leftLine = new TLine(dmSigmas["left"], yMin, dmSigmas["left"], yMaxCanvas);
    rightLine = new TLine(dmSigmas["right"], yMin, dmSigmas["right"], yMaxCanvas);
    leftLine->SetLineColor(kBlue);
    leftLine->SetLineStyle(2);
    leftLine->SetLineWidth(2);
    rightLine->SetLineColor(kBlue);
    rightLine->SetLineStyle(2);
    rightLine->SetLineWidth(2);

    leftLine->Draw("same");
    rightLine->Draw("same");
    latex.SetTextSize(0.035);
    latex.SetTextFont(42);
    latex.SetTextColor(kBlack);
    left_sigma = dmSigmas["bestfit"] - dmSigmas["left"];
    right_sigma = dmSigmas["right"] - dmSigmas["bestfit"];
    latex.DrawLatex(dmSigmas["left"] - 0.000029, yMaxCanvas * 0.9, Form("#Deltam^{2}_{21} = (%.3f^{+%.3f}_{-%.3f}) #times10^{-5} eV^{2}", dmSigmas["bestfit"] * 1E5, right_sigma * 1E5, left_sigma * 1E5));

    pathObj.replace_filename("deltam21LLHDiff.pdf");
    c3->SaveAs(pathObj.string().c_str());
    outfile->cd();
    c3->Write("c3");
    Profile.second->Write("hDeltam");

    outfile->WriteObject(&dmSigmas, "dmSigmas");
    outfile->WriteObject(&thSigmas, "thSigmas");
}
