#include <TFile.h>
#include <TTree.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TStyle.h>

/* ///////////////////////////////////////////////////////////////////
///
/// Script for plotting the postfit energy distributions after fixed
/// oscillation fits.
/// First open the root file (made by makeFixedOscTree) and find the
/// minimum LLH entry.
/// Then go into the corresponding fit directory and load up the
/// histograms from scaled_dists
/// Loop through these and plot, and also add to groups for a separate
/// groups plot.
/// Canvases are saved to .pdf and .root
///
/////////////////////////////////////////////////////////////////// */

void plotFixedOscDist(const char *filename = "fit_results.root")
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
        std::cerr << "Error: Could not find TTree 'fitResults' in file " << filename << std::endl;
        file->Close();
        return;
    }

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
            bestTheta = branchValues["theta12"];
            bestDeltaM = branchValues["deltam21"];
        }
    }

    // Declare group histos and stacks
    TH1D *hReactor;
    TH1D *hGeo;
    TH1D *hAlpha;
    TH1D *hOther;
    TH1D *hTotal;
    TH1D *hData;
    THStack *hStack = new THStack("hStack", "Stacked;Reconstructed Energy, MeV;");
    THStack *hGroupStack = new THStack("hGroupStack", "GroupStacked;Reconstructed Energy, MeV;");

    // Map of filenames to plot labels
    std::map<std::string, std::string> labelMap;
    labelMap["alphan_CScatter.root"] = "#alpha C Scatter ";
    labelMap["alphan_OExcited.root"] = "#alpha O excited";
    labelMap["alphan_PRecoil.root"] = "#alpha P Recoil";
    labelMap["data.root"] = "Data";
    labelMap["geonu_Th.root"] = "Geo Th";
    labelMap["geonu_U.root"] = "Geo U";
    labelMap["postfitdist.root"] = "Total MC with Systs.";
    labelMap["reactor_nubar.root"] = "Reactor #bar{#nu} ";
    labelMap["sideband.root"] = "Sideband";

    // Directory of best LLH fit
    std::filesystem::path dirPath(filename);
    dirPath.replace_filename("");
    std::ostringstream directory;
    directory << dirPath.string() << "/th" << std::fixed << std::setprecision(2) << bestTheta
              << "/th" << std::fixed << std::setprecision(2) << bestTheta
              << "_dm" << std::fixed << std::setprecision(8) << bestDeltaM
              << "/scaled_dists/";

    // Declare canvas and legends
    TCanvas *c1 = new TCanvas("c1", "Distributions", 800, 600);
    gPad->SetGrid();
    gStyle->SetOptStat(0);

    TLegend *t1 = new TLegend(0.5, 0.5, 0.85, 0.85);
    t1->SetLineWidth(2);
    TLegend *t2 = new TLegend(0.5, 0.5, 0.85, 0.85);
    t2->SetLineWidth(2);

    // Define colours for histograms
    std::vector<int> lineColors = {kRed + 3, kRed + 2, kRed + 1, kMagenta + 2, kMagenta + 4, kBlue, kGreen + 3, kGray, kBlack};
    std::vector<int> fillColors = {kRed - 1, kRed - 2, kRed - 9, kMagenta - 8, kMagenta - 5, kBlue - 9, kGreen - 5, kGray, kBlack};
    int colorIndex = 0;

    bool firstHist = true;
    bool firstDraw = true;

    // Go into scaled dists and loop over each file
    std::cout << directory.str() << std::endl;
    for (const auto &entry : std::filesystem::directory_iterator(directory.str()))
    {

        std::string filePath = entry.path().string();
        std::filesystem::path pathObj(filePath);

        // Check if it's a .root file
        if (filePath.find(".root") == std::string::npos)
            continue;

        TFile *file = TFile::Open(filePath.c_str(), "READ");
        if (!file || file->IsZombie())
        {
            std::cerr << "Error opening file: " << filePath << std::endl;
            delete file;
            continue;
        }

        // Get histogram
        TH1D *h1 = nullptr;
        file->GetObject("oxsx_saved", h1);

        if (!h1)
        {
            std::cerr << "No histogram 'oxsx_saved' found in file: " << filePath << std::endl;
            file->Close();
            delete file;
            continue;
        }

        // Set histogram properties
        h1->SetLineColor(lineColors[colorIndex % lineColors.size()]);
        h1->SetLineWidth(2);
        h1->GetYaxis()->SetRangeUser(0, 1.6);

        // Initialise group histograms
        if (firstHist)
        {
            hData = (TH1D *)h1->Clone("hData");
            hData->Reset();
            hTotal = (TH1D *)h1->Clone("hTotal");
            hTotal->Reset();
            hReactor = (TH1D *)h1->Clone("hReactor");
            hReactor->Reset();
            hGeo = (TH1D *)h1->Clone("hGeo");
            hGeo->Reset();
            hAlpha = (TH1D *)h1->Clone("hAlpha");
            hAlpha->Reset();
            hOther = (TH1D *)h1->Clone("hOther");
            hOther->Reset();
            firstHist = false;
        }

        // If we're data or summed MC then do not add to the stack
        if (pathObj.filename() != "data.root" && pathObj.filename() != "postfitdist.root")
        {
            h1->SetFillColor(fillColors[colorIndex % fillColors.size()]);
            hStack->Add(h1);
            t1->AddEntry(h1, labelMap[pathObj.filename()].c_str(), "f");
        }
        else
        {
            t1->AddEntry(h1, labelMap[pathObj.filename()].c_str(), "l");
            if (firstDraw)
            {
                h1->Draw("HIST");
                h1->GetYaxis()->SetRangeUser(0, 1.6);
                firstDraw = 0;
            }
            else
                h1->Draw("same HIST");
        }

        colorIndex++;

        // Now add distributions to groups
        if (pathObj.filename().string() == "data.root")
        {
            hData->Add(h1);
            hData->SetLineColor(kBlack);
        }
        else if (pathObj.filename().string() == "postfitdist.root")
        {
            hTotal->Add(h1);
            hTotal->SetLineColor(kGray + 2);
        }
        else if (pathObj.filename().string() == "reactor_nubar.root")
        {
            hReactor->Add(h1);
            hReactor->SetLineColor(kBlue);
            hReactor->SetFillColor(kBlue - 9);
        }
        else if (pathObj.filename().string() == "geonu_U.root" || pathObj.filename().string() == "geonu_Th.root")
        {
            hGeo->Add(h1);
            hGeo->SetLineColor(kMagenta + 2);
            hGeo->SetFillColor(kMagenta - 8);
        }
        else if (pathObj.filename().string() == "alphan_CScatter.root" || pathObj.filename().string() == "alphan_OExcited.root" || pathObj.filename().string() == "alphan_PRecoil.root")
        {
            hAlpha->Add(h1);
            hAlpha->SetLineColor(kRed + 1);
            hAlpha->SetFillColor(kRed - 9);
        }
        else
        {
            hOther->Add(h1);
            hOther->SetLineColor(kGreen + 3);
            hOther->SetFillColor(kGreen - 5);
        }
    }

    // Draw the legend
    t1->SetLineWidth(2);
    hStack->Draw("same");
    t1->Draw();

    // Update the canvas
    c1->Update();
    std::filesystem::path pathObj(filename);
    pathObj.replace_filename("dist.pdf");
    c1->SaveAs(pathObj.string().c_str());
    pathObj.replace_filename("dist.root");
    c1->SaveAs(pathObj.string().c_str());

    // Add each group histos to stack
    TCanvas *c2 = new TCanvas("c2", "Groups", 800, 600);
    gPad->SetGrid();
    gStyle->SetOptStat(0);

    t2->AddEntry(hData, "Data", "l");
    t2->AddEntry(hTotal, "Total MC with Systs.", "l");
    hGroupStack->Add(hReactor);
    t2->AddEntry(hReactor, "Reactor #bar{#nu}", "f");
    hGroupStack->Add(hAlpha);
    t2->AddEntry(hAlpha, "#alpha-n", "f");
    hGroupStack->Add(hGeo);
    t2->AddEntry(hGeo, "Geo #bar{#nu}", "f");
    hGroupStack->Add(hOther);
    t2->AddEntry(hOther, "Other", "f");

    hGroupStack->Draw();
    hData->Draw("same");
    hTotal->Draw("same");
    t2->Draw();

    // Update the canvas
    c2->Update();
    pathObj.replace_filename("distGroup.pdf");
    c2->SaveAs(pathObj.string().c_str());
    pathObj.replace_filename("distGroup.root");
    c2->SaveAs(pathObj.string().c_str());
}
