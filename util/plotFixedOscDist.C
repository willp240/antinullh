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
///
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
    TH1D *hData;
    THStack *hStack = new THStack("hStack", ";Reconstructed Energy, MeV;Events");
    THStack *hGroupStack = new THStack("hGroupStack", ";Reconstructed Energy, MeV;Events");
    std::map<std::string, TH1D *> histMap;

    // Make map of param names to plot labels
    std::map<std::string, std::string> labelMap;
    std::vector<std::string> *paramNames = nullptr;
    std::vector<std::string> *labelsVec = nullptr;
    file->GetObject("param_names", paramNames);
    file->GetObject("tex_labels", labelsVec);

    for (int iParam = 0; iParam < paramNames->size(); iParam++)
    {
        labelMap[paramNames->at(iParam)] = labelsVec->at(iParam);
    }
    labelMap["data"] = "Data";

    std::vector<std::string> paramOrder{"data", "reactor_nubar", "alphan_PRecoil", "alphan_CScatter",
                                       "alphan_OExcited", "geonu_Th", "geonu_U", "sideband"};

    // Directory of best LLH fit
    std::filesystem::path dirPath(filename);
    dirPath.replace_filename("");
    std::ostringstream directory;
    directory << dirPath.string() << "/th" << std::fixed << std::setprecision(2) << bestTheta
              << "/th" << std::fixed << std::setprecision(2) << bestTheta
              << "_dm" << std::fixed << std::setprecision(8) << bestDeltaM
              << "/scaled_dists/";

    TLegend *t1 = new TLegend(0.5, 0.48, 0.85, 0.85);
    t1->SetLineWidth(2);
    TLegend *t2 = new TLegend(0.5, 0.5, 0.85, 0.85);
    t2->SetLineWidth(2);

    // Define colours for histograms
    std::vector<int> lineColours = {kBlack, kBlue, kMagenta + 2, kMagenta + 4, kRed + 3, kRed + 2, kRed + 1, kGreen + 3};
    std::vector<int> fillColours = {kBlack, kBlue - 9, kMagenta - 8, kMagenta - 5, kRed - 1, kRed - 2, kRed - 9, kGreen - 5};
    int colourIndex = 0;

    // Go into scaled dists and loop over each file
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
        pathObj.replace_extension("");
        histMap[pathObj.filename()] = h1;
    }

    // Intialise group histos
    hData = (TH1D *)histMap["data"]->Clone("hData");
    hData->Reset();
    hReactor = (TH1D *)histMap["data"]->Clone("hReactor");
    hReactor->Reset();
    hGeo = (TH1D *)histMap["data"]->Clone("hGeo");
    hGeo->Reset();
    hAlpha = (TH1D *)histMap["data"]->Clone("hAlpha");
    hAlpha->Reset();
    hOther = (TH1D *)histMap["data"]->Clone("hOther");
    hOther->Reset();

    hData->SetLineColor(kBlack);
    hData->SetLineWidth(2);
    hReactor->SetLineColor(kBlue);
    hReactor->SetFillColor(kBlue - 9);
    hReactor->SetLineWidth(2);
    hGeo->SetLineColor(kMagenta + 2);
    hGeo->SetFillColor(kMagenta - 8);
    hGeo->SetLineWidth(2);
    hAlpha->SetLineColor(kRed + 1);
    hAlpha->SetFillColor(kRed - 9);
    hAlpha->SetLineWidth(2);
    hOther->SetLineColor(kGreen + 3);
    hOther->SetFillColor(kGreen - 5);
    hOther->SetLineWidth(2);

    colourIndex = 0;

    // Now loop over files, getting histos, and adding to the appropriate groups
    for (int iFile = 0; iFile < paramOrder.size(); iFile++)
    {
        histMap[paramOrder.at(iFile)]->SetLineColor(lineColours[colourIndex % lineColours.size()]);
        histMap[paramOrder.at(iFile)]->SetLineWidth(2);
        histMap[paramOrder.at(iFile)]->GetYaxis()->SetRangeUser(0, 1.6);

        if (paramOrder.at(iFile) == "data")
        {
            hData->Add(histMap[paramOrder.at(iFile)]);
            t1->AddEntry(histMap[paramOrder.at(iFile)], labelMap[paramOrder.at(iFile)].c_str(), "l");
        }
        else if (paramOrder.at(iFile) == "reactor_nubar")
        {
            hReactor->Add(histMap[paramOrder.at(iFile)]);
            hStack->Add(histMap[paramOrder.at(iFile)]);
            histMap[paramOrder.at(iFile)]->SetFillColor(fillColours[colourIndex % fillColours.size()]);
            t1->AddEntry(histMap[paramOrder.at(iFile)], labelMap[paramOrder.at(iFile)].c_str(), "f");
            histMap.erase(paramOrder.at(iFile));
        }
        else if (paramOrder.at(iFile) == "geonu_U" || paramOrder.at(iFile) == "geonu_Th")
        {
            hGeo->Add(histMap[paramOrder.at(iFile)]);
            hStack->Add(histMap[paramOrder.at(iFile)]);
            histMap[paramOrder.at(iFile)]->SetFillColor(fillColours[colourIndex % fillColours.size()]);
            t1->AddEntry(histMap[paramOrder.at(iFile)], labelMap[paramOrder.at(iFile)].c_str(), "f");
            histMap.erase(paramOrder.at(iFile));
        }
        else if (paramOrder.at(iFile) == "alphan_CScatter" || paramOrder.at(iFile) == "alphan_OExcited" || paramOrder.at(iFile) == "alphan_PRecoil")
        {
            hAlpha->Add(histMap[paramOrder.at(iFile)]);
            hStack->Add(histMap[paramOrder.at(iFile)]);
            histMap[paramOrder.at(iFile)]->SetFillColor(fillColours[colourIndex % fillColours.size()]);
            t1->AddEntry(histMap[paramOrder.at(iFile)], labelMap[paramOrder.at(iFile)].c_str(), "f");
            histMap.erase(paramOrder.at(iFile));
        }
        else
        {
            hOther->Add(histMap[paramOrder.at(iFile)]);
            hStack->Add(histMap[paramOrder.at(iFile)]);
            histMap[paramOrder.at(iFile)]->SetFillColor(fillColours[colourIndex % fillColours.size()]);
            t1->AddEntry(histMap[paramOrder.at(iFile)], labelMap[paramOrder.at(iFile)].c_str(), "f");
            histMap.erase(paramOrder.at(iFile));
        }

        colourIndex++;
    }

    colourIndex = 0;
    // Loop over what's left in histMap in case there were any files not defined in labelMap
    for (std::map<std::string, TH1D *>::iterator it = histMap.begin(); it != histMap.end(); ++it)
    {
        if (it->first == "data" || it->first == "postfitdist")
            continue;

        hOther->Add(histMap[it->first]);
        histMap[it->first]->SetFillColor(colourIndex + 30);
        t1->AddEntry(histMap[it->first], it->first.c_str(), "f");
        colourIndex++;
    }

    // Draw stack of all event types
    TCanvas *c1 = new TCanvas("c1", "Stacked", 1000, 600);
    gPad->SetGrid();
    gStyle->SetOptStat(0);
    c1->SetFrameLineWidth(2);

    hStack->SetMaximum(1.6);
    hStack->Draw("");
    histMap["data"]->Draw("histsame");
    hStack->GetXaxis()->SetTitleFont(42);
    hStack->GetYaxis()->SetTitleFont(42);
    hStack->GetXaxis()->SetLabelFont(42);
    hStack->GetYaxis()->SetLabelFont(42);
    t1->SetTextFont(42);
    t1->Draw();
    c1->Update();

    std::filesystem::path pathObj(filename);
    pathObj.replace_filename("dist.pdf");
    c1->SaveAs(pathObj.string().c_str());
    pathObj.replace_filename("dist.root");
    c1->SaveAs(pathObj.string().c_str());

    // Draw stack of grouped event types
    TCanvas *c2 = new TCanvas("c2", "Groups", 1000, 600);
    gPad->SetGrid();
    gStyle->SetOptStat(0);
    c2->SetFrameLineWidth(2);

    t2->AddEntry(hData, "Data", "l");
    hGroupStack->Add(hReactor);
    t2->AddEntry(hReactor, "Reactor #bar{#nu}", "f");
    hGroupStack->Add(hAlpha);
    t2->AddEntry(hAlpha, "#alpha-n", "f");
    hGroupStack->Add(hGeo);
    t2->AddEntry(hGeo, "Geo #bar{#nu}", "f");
    hGroupStack->Add(hOther);
    t2->AddEntry(hOther, "Other", "f");

    hGroupStack->SetMaximum(1.6);
    hGroupStack->Draw("");
    histMap["data"]->Draw("histsame");
    hGroupStack->GetXaxis()->SetTitleFont(42);
    hGroupStack->GetYaxis()->SetTitleFont(42);
    hGroupStack->GetXaxis()->SetLabelFont(42);
    hGroupStack->GetYaxis()->SetLabelFont(42);
    t2->SetTextFont(42);
    t2->Draw();
    c2->Update();

    pathObj.replace_filename("distGroup.pdf");
    c2->SaveAs(pathObj.string().c_str());
    pathObj.replace_filename("distGroup.root");
    c2->SaveAs(pathObj.string().c_str());
}
