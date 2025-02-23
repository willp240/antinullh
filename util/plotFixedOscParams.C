#include <TFile.h>
#include <TTree.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TStyle.h>

/* ///////////////////////////////////////////////////////////////////
///
/// Script for plotting post fit parameter values and prefit
/// constraints relative to nominal values for fixed oscillation fits.
///
/// The user inputs the root file made by makeFixedOscTree, and it
/// first finds the entry with the lowest LLH. All the branches for
/// this entry are read into a map, and the nominal values and
/// constraints vectors are read in from the input file too. Those
/// vectors are sorted into the order we want to plot (so this will
/// need to be updated when parameters change).
///
/// The plot is drawn and the canvas is saved as a root file and
/// pdf.
///
/////////////////////////////////////////////////////////////////// */

/// Function to sort vectors of names, nominals, and constraints into the order we want to plot them in.
/// Also made is a vector of Latex parameter labels (we could (should?) read these from a config)
void sortVectors(std::vector<std::string> *&namesVec, std::vector<std::string> *&labelsVec, std::vector<double> *&nomsVec, std::vector<double> *&constrMeansVec, std::vector<double> *&constrErrsVec)
{

    std::vector<double> *tempNomsVec = new std::vector<double>(namesVec->size(), 0.0);
    std::vector<double> *tempConstrMeansVec = new std::vector<double>(namesVec->size(), 0.0);
    std::vector<double> *tempConstrErrsVec = new std::vector<double>(namesVec->size(), 0.0);
    std::vector<std::string> *tempLabelsVec = new std::vector<std::string>(namesVec->size(), "");

    // This is the order we want to plot them in (osc, signal, geo, alpha n, other bgs, systematics)
    std::vector<std::string> *tempNamesVec = new std::vector<std::string>{"deltam21",
                                                                          "theta12",
                                                                          "reactor_nubar",
                                                                          "geonu_U",
                                                                          "geonu_Th",
                                                                          "alphan_CScatter",
                                                                          "alphan_OExcited",
                                                                          "alphan_PRecoil",
                                                                          "sideband",
                                                                          "energy_scale",
                                                                          "energy_conv",
                                                                          "birks_constant",
                                                                          "p_recoil_energy_scale"};

    // Now loop over the new vector (that's already in the correct order), and find the index of the same name
    // in the old vector, and set the noms and constraints to be the ones for that index
    for (int iTempName = 0; iTempName < tempNamesVec->size(); iTempName++)
    {
        for (int iName = 0; iName < namesVec->size(); iName++)
        {
            if (tempNamesVec->at(iTempName) == namesVec->at(iName))
            {
                tempNomsVec->at(iTempName) = nomsVec->at(iName);
                tempConstrMeansVec->at(iTempName) = constrMeansVec->at(iName);
                tempConstrErrsVec->at(iTempName) = constrErrsVec->at(iName);
                tempLabelsVec->at(iTempName) = labelsVec->at(iName);
                break;
            }
        }
    }

    namesVec = tempNamesVec;
    labelsVec = tempLabelsVec;
    nomsVec = tempNomsVec;
    constrMeansVec = tempConstrMeansVec;
    constrErrsVec = tempConstrErrsVec;
}

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
    int bestLLHEntry = 999;

    // Loop over entries in the tree to find the lowest LLH step
    Long64_t nEntries = tree->GetEntries();
    for (Long64_t i = 0; i < nEntries; ++i)
    {
        tree->GetEntry(i);
        if (branchValues["LLH"] < minLLH)
        {
            minLLH = branchValues["LLH"];
            bestLLHEntry = i;
        }
    }

    tree->GetEntry(bestLLHEntry);

    std::cout << "Best LLH: " << branchValues["LLH"] << std::endl;
    std::cout << "Best Entry: " << bestLLHEntry << std::endl;
    std::cout << "Best theta: " << branchValues["theta12"] << std::endl;
    std::cout << "Best deltam: " << branchValues["deltam21"] << std::endl;

    // Retrieve the vectors from the file: nominal, constr means, constr err, constr names
    std::vector<std::string> *paramNames = nullptr;
    std::vector<double> *nomVals = nullptr;
    std::vector<double> *constrMeans = nullptr;
    std::vector<double> *constrErr = nullptr;
    std::vector<std::string> *labelsVec = nullptr;

    file->GetObject("param_names", paramNames);
    file->GetObject("param_asimov_values", nomVals);
    file->GetObject("constr_mean_values", constrMeans);
    file->GetObject("constr_sigma_values", constrErr);
    file->GetObject("tex_labels", labelsVec);

    // Reorder vectors to the order we want to plot them
    sortVectors(paramNames, labelsVec, nomVals, constrMeans, constrErr);

    TH1D *hNom = new TH1D("hNominal", "Relative Nominal Values", paramNames->size(), 0, paramNames->size() - 1);
    TH1D *hConstr = new TH1D("hConstr", "Constraints Relative to Nominals", paramNames->size(), 0, paramNames->size() - 1);
    TH1D *hPostFit = new TH1D("hPostFit", "Postfit Values Relative to Nominals", paramNames->size(), 0, paramNames->size() - 1);

    for (int iParam = 0; iParam < paramNames->size(); iParam++)
    {
        // For reactor nubar, rate in fit config is unoscillated, but postfit value is oscillated.
        // So we need to multiply by the ratio saved in outputted fit file
        if (paramNames->at(iParam) == "reactor_nubar")
        {
            nomVals->at(iParam) = nomVals->at(iParam) * branchValues["reactor_ratio"];
            constrMeans->at(iParam) = constrMeans->at(iParam) * branchValues["reactor_ratio"];
            constrErr->at(iParam) = constrErr->at(iParam) * branchValues["reactor_ratio"];
        }

        hNom->SetBinContent(iParam + 1, nomVals->at(iParam) / nomVals->at(iParam));
        hConstr->GetXaxis()->SetBinLabel(iParam + 1, labelsVec->at(iParam).c_str());
        hConstr->SetBinContent(iParam + 1, constrMeans->at(iParam) / nomVals->at(iParam));
        hConstr->SetBinError(iParam + 1, constrErr->at(iParam) / nomVals->at(iParam));
        hPostFit->SetBinContent(iParam + 1, branchValues[paramNames->at(iParam)] / nomVals->at(iParam));
        // hPostFit->SetBinError(iParam + 1, branchValues[paramNames->at(iParam) + "_err"] / nomVals->at(iParam));
        // TODO: If fit valid is 0, we get 0 bars and the plot goes wacky. Just set to a value now to get a marker, but we should fix this
        hPostFit->SetBinError(iParam + 1, 0.1);
    }

    // Draw the histograms
    TCanvas *c1 = new TCanvas("c1", "Params", 800, 600);
    c1->SetBottomMargin(0.23);
    gPad->SetFrameLineWidth(2);
    gStyle->SetOptStat(0);
    gPad->SetGrid(1);

    hNom->SetLineColor(kGreen);
    hNom->SetLineWidth(2);
    hConstr->SetLineColor(kBlue);
    hConstr->SetLineWidth(2);
    hPostFit->SetLineColor(kRed);
    hPostFit->SetLineWidth(2);
    // TODO: We'll want these when we have proper error bars
    // hPostFit->SetMarkerStyle(2);
    // hPostFit->SetMarkerSize(4);
    // hPostFit->SetMarkerColor(kRed);

    hConstr->GetYaxis()->SetRangeUser(0, 2);
    hConstr->GetYaxis()->SetTitle("Relative to Nominal");
    //hConstr->GetXaxis()->SetLabelOffset(0.007);
    hConstr->GetXaxis()->CenterLabels(true);
    hConstr->SetTitle("");
    hConstr->LabelsOption("v");

    hConstr->GetXaxis()->SetTitleFont(42);
    hConstr->GetYaxis()->SetTitleFont(42);
    hConstr->GetXaxis()->SetLabelFont(42);
    hConstr->GetYaxis()->SetLabelFont(42);
    hConstr->SetTitleFont(42);

    hConstr->Draw("E1");
    hPostFit->Draw("E1same");

    TLegend *t1 = new TLegend(0.65, 0.73, 0.85, 0.88);
    t1->AddEntry(hConstr, "Prefit", "l");
    t1->AddEntry(hPostFit, "Postfit", "l");
    t1->SetLineWidth(2);
    t1->SetTextFont(42);
    t1->Draw();

    // Save plot as image and rootfile
    std::filesystem::path pathObj(filename);
    pathObj.replace_filename("params.pdf");
    c1->SaveAs(pathObj.string().c_str());
    pathObj.replace_filename("params.root");
    c1->SaveAs(pathObj.string().c_str());
}
