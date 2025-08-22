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
/// Script for plotting the postfit energy distributions of three
/// fixed oscillation fits.
/// First open the root files (made by makeFixedOscTree) and find the
/// minimum LLH entries.
/// Then go into the corresponding fit directories and load up the
/// histograms from scaled_dists
///
/// Plot the data from the first fit, and the total fitted distribution
/// from each fit
///
/// Canvases are saved to .pdf and .root
///
/////////////////////////////////////////////////////////////////// */

std::ostringstream GetBestFitDir(const char *filename, TTree *tree)
{

    // Get the list of branches
    TObjArray *branches = tree->GetListOfBranches();
    std::map<std::string, double> branchValues;

    // Create a map of branch names to variables
    std::map<std::string, double *> branchPointers;
    for (int i = 0; i < branches->GetEntries(); ++i)
    {
        std::string branchName = branches->At(i)->GetName();
        branchValues[branchName] = 0.0;
        branchPointers[branchName] = &branchValues[branchName];
        tree->SetBranchAddress(branchName.c_str(), branchPointers[branchName]);
    }

    // Check the form of theta 12
    std::string theta12name;
    if (tree->GetBranch("theta12"))
        theta12name = "theta12";
    else if (tree->GetBranch("sintheta12"))
        theta12name = "sintheta12";
    else if (tree->GetBranch("sinsqtheta12"))
        theta12name = "sinsqtheta12";

    double minLLH = 1e9;
    double bestTheta = 0;
    double bestDeltaM = 0;

    // Loop over entries in the tree
    Long64_t nEntries = tree->GetEntries();
    for (Long64_t i = 0; i < nEntries; ++i)
    {
        tree->GetEntry(i);
        if (branchValues["LLH"] < minLLH)
        {
            minLLH = branchValues["LLH"];
            bestTheta = branchValues[theta12name];
            bestDeltaM = branchValues["deltam21"];
        }
    }

    // Directory of best LLH fit
    std::filesystem::path dirPath(filename);
    dirPath.replace_filename("");
    std::ostringstream directory;
    directory << dirPath.string() << "/th" << std::fixed << std::setprecision(3) << bestTheta
              << "/th" << std::fixed << std::setprecision(3) << bestTheta
              << "_dm" << std::fixed << std::setprecision(8) << bestDeltaM
              << "/postfit_dists/";

    return directory;
}

void compare3FixedOscDists(const char *filename1, const char *filename2, const char *filename3, std::string label1, std::string label2, std::string label3)
{

    // Open the ROOT files
    TFile *file1 = TFile::Open(filename1, "READ");
    if (!file1 || file1->IsZombie())
    {
        std::cerr << "Error: Could not open file " << filename1 << std::endl;
        return;
    }
    TFile *file2 = TFile::Open(filename2, "READ");
    if (!file2 || file2->IsZombie())
    {
        std::cerr << "Error: Could not open file " << filename2 << std::endl;
        return;
    }
    TFile *file3 = TFile::Open(filename3, "READ");
    if (!file3 || file3->IsZombie())
    {
        std::cerr << "Error: Could not open file " << filename3 << std::endl;
        return;
    }

    // Get the TTrees
    TTree *tree1 = nullptr;
    file1->GetObject("fitResults", tree1);
    if (!tree1)
    {
        std::cerr << "Error: Could not find TTree 'fitResults' in file " << filename1 << std::endl;
        file1->Close();
        return;
    }
    TTree *tree2 = nullptr;
    file2->GetObject("fitResults", tree2);
    if (!tree2)
    {
        std::cerr << "Error: Could not find TTree 'fitResults' in file " << filename2 << std::endl;
        file2->Close();
        return;
    }
    TTree *tree3 = nullptr;
    file3->GetObject("fitResults", tree3);
    if (!tree3)
    {
        std::cerr << "Error: Could not find TTree 'fitResults' in file " << filename3 << std::endl;
        file3->Close();
        return;
    }

    std::ostringstream bestfitdir1 = GetBestFitDir(filename1, tree1);
    std::ostringstream bestfitdir2 = GetBestFitDir(filename2, tree2);
    std::ostringstream bestfitdir3 = GetBestFitDir(filename3, tree3);

    std::filesystem::path datapath = std::filesystem::path(bestfitdir1.str() + "/../data.root");
    TFile *datafile1 = TFile::Open(datapath.c_str(), "READ");
    TH1D *hdata = nullptr;
    datafile1->GetObject("oxsx_saved", hdata);

    TFile *scaledfile1 = TFile::Open((bestfitdir1.str() + "/postfitdist.root").c_str(), "READ");
    TH1D *hscaled1 = nullptr;
    scaledfile1->GetObject("oxsx_saved", hscaled1);

    TFile *scaledfile2 = TFile::Open((bestfitdir2.str() + "/postfitdist.root").c_str(), "READ");
    TH1D *hscaled2 = nullptr;
    scaledfile2->GetObject("oxsx_saved", hscaled2);

    TFile *scaledfile3 = TFile::Open((bestfitdir3.str() + "/postfitdist.root").c_str(), "READ");
    TH1D *hscaled3 = nullptr;
    scaledfile3->GetObject("oxsx_saved", hscaled3);

    TLegend *t1 = new TLegend(0.5, 0.48, 0.85, 0.85);
    t1->SetLineWidth(2);
    t1->AddEntry(hdata, "Data", "l");
    t1->AddEntry(hscaled1, label1.c_str(), "l");
    t1->AddEntry(hscaled2, label2.c_str(), "l");
    t1->AddEntry(hscaled3, label3.c_str(), "l");

    hdata->SetLineColor(kBlack);
    hdata->SetLineWidth(4);
    hscaled1->SetLineColor(kBlue + 2);
    hscaled1->SetLineWidth(4);
    hscaled1->SetLineStyle(2);
    hscaled2->SetLineColor(kRed + 2);
    hscaled2->SetLineWidth(4);
    hscaled2->SetLineStyle(3);
    hscaled3->SetLineColor(kGreen + 2);
    hscaled3->SetLineWidth(4);
    hscaled3->SetLineStyle(4);

    // Draw stack of all event types
    TCanvas *c1 = new TCanvas("c1", "c1", 1000, 600);
    gPad->SetGrid();
    gStyle->SetOptStat(0);
    c1->SetFrameLineWidth(2);
    double uppermin = 0.3;
    TPad *lower = new TPad("lower", "pad", 0, 0, 1, uppermin);
    TPad *upper = new TPad("upper", "pad", 0, uppermin, 1, 1);
    upper->SetBottomMargin(0.04);
    lower->SetTopMargin(0.06);
    lower->SetBottomMargin(0.4);
    upper->SetFrameLineWidth(2);
    lower->SetFrameLineWidth(2);
    upper->SetGrid();
    lower->SetGrid();
    upper->Draw();
    lower->Draw();
    c1->cd();

    upper->cd();
    hdata->SetMaximum(1.6);
    hdata->Draw("");
    hscaled1->Draw("histsame");
    hscaled2->Draw("histsame");
    hscaled3->Draw("histsame");
    hdata->GetXaxis()->SetLabelOffset(1.2);
    hdata->GetXaxis()->SetTitleFont(42);
    hdata->GetYaxis()->SetTitleFont(42);
    hdata->GetXaxis()->SetLabelFont(42);
    hdata->GetYaxis()->SetLabelFont(42);
    hdata->GetXaxis()->SetLabelSize(0.06);
    hdata->GetYaxis()->SetLabelSize(0.06);
    hdata->GetXaxis()->SetTitleSize(0.06);
    hdata->GetYaxis()->SetTitleSize(0.06);
    hdata->GetYaxis()->SetTitleOffset(0.8);
    t1->SetTextFont(42);
    t1->Draw();
    c1->Update();

    TH1D *hscaled1Div = (TH1D *)hscaled1->Clone("hscaled1Div");
    TH1D *hscaled2Div = (TH1D *)hscaled2->Clone("hscaled2Div");
    TH1D *hscaled3Div = (TH1D *)hscaled3->Clone("hscaled3Div");

    lower->cd();
    hscaled1Div->Divide(hdata);
    hscaled2Div->Divide(hdata);
    hscaled3Div->Divide(hdata);
    hscaled1Div->GetYaxis()->SetRangeUser(0.9, 1.1);
    hscaled1Div->GetXaxis()->SetTitleFont(42);
    hscaled1Div->GetYaxis()->SetTitleFont(42);
    hscaled1Div->GetXaxis()->SetLabelFont(42);
    hscaled1Div->GetYaxis()->SetLabelFont(42);
    hscaled1Div->GetXaxis()->SetLabelSize(0.14);
    hscaled1Div->GetYaxis()->SetLabelSize(0.14);
    hscaled1Div->GetXaxis()->SetTitleSize(0.14);
    hscaled1Div->GetYaxis()->SetTitleSize(0.14);
    hscaled1Div->GetYaxis()->SetTitle("MC/Data");
    hscaled1Div->GetYaxis()->SetTitleOffset(0.33);
    hscaled1Div->GetYaxis()->SetNdivisions(404);
    hscaled1Div->GetYaxis()->SetTickLength(0.05);
    hscaled1Div->GetXaxis()->SetTickLength(0.07);
    hscaled1Div->GetYaxis()->ChangeLabel(2, -1, -1, -1, -1, -1, " ");
    hscaled1Div->GetYaxis()->ChangeLabel(4, -1, -1, -1, -1, -1, " ");
    hscaled1Div->SetLineWidth(2);
    hscaled1Div->Draw();
    hscaled2Div->Draw("same");
    hscaled3Div->Draw("same");
    c1->Update();

    struct stat st = {0};
    std::filesystem::path pathObj(filename1);
    pathObj.replace_filename("plots/");
    if (stat(pathObj.string().c_str(), &st) == -1)
        mkdir(pathObj.string().c_str(), 0700);
    pathObj.replace_filename("distcomp.pdf");
    c1->SaveAs(pathObj.string().c_str());
    pathObj.replace_filename("distcomp.root");
    c1->SaveAs(pathObj.string().c_str());
}
