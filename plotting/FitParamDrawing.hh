#pragma once

#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TColor.h>
#include <TMatrixT.h>
#include <TGraphAsymmErrors.h>
#include <TNamed.h>
#include <algorithm>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <map>
#include <string>
#include <vector>
#include <unordered_map>
#include <sys/stat.h>

namespace antinupostfit
{
// Hard code reactor parameter names to make sure we use the right ratio for each
// Inelegant, but we hard code the parameter order anyway
const std::string reactorpar1 = "reactor_nubar";
const std::string reactorpar2 = "reactor_nubar2";

/// Function to sort vectors of names, nominals, and constraints into the order we want to plot them in
inline void sortVectors(std::vector<std::string> *&namesVec, std::vector<double> *&errVec, std::vector<std::string> *&labelsVec, std::vector<double> *&nomsVec,
                 std::vector<double> *&constrMeansVec, std::vector<double> *&constrErrsVec, TMatrixT<double> *&covMatrix, std::string theta12name, bool corr, bool alphanclass)
{

    int nParams = namesVec->size();
    std::vector<std::string> *orderedNamesVec;

    if (alphanclass)
    {
        if (!corr)
        {
            orderedNamesVec = new std::vector<std::string>{
                "deltam21",
                theta12name,
                "reactor_nubar",
                "geonu_U",
                "geonu_Th",
                "alphan_CScatter",
                "alphan_OExcited",
                "alphan_PRecoil",
                "bipolike",
                "atmospheric",
                "energy_scale",
                "energy_conv",
                "birks_constant",
                "p_recoil_energy_scale",
                "reactor_nubar2",
                "geonu_U2",
                "geonu_Th2",
                "alphan_CScatter2",
                "alphan_OExcited2",
                "alphan_PRecoil2",
                "bipolike2",
                "atmospheric2",
                "energy_scale2",
                "energy_conv2",
                "birks_constant2",
                "p_recoil_energy_scale2",
                "class_a_ppo",
                "class_a_bismsb",
                "class_s_ppo",
                "class_s_bismsb"};
        }
        else
        {
            orderedNamesVec = new std::vector<std::string>{
                "deltam21",
                theta12name,
                "reactor_nubar_norm",
                "geonu_U_norm",
                "geonu_Th_norm",
                "alphan_CScatter_norm",
                "alphan_OExcited_norm",
                "alphan_PRecoil_norm",
                "bipolike_norm",
                "atmospheric_norm",
                "energy_scale",
                "energy_conv",
                "birks_constant",
                "p_recoil_energy_scale",
                "energy_scale2",
                "energy_conv2",
                "birks_constant2",
                "p_recoil_energy_scale2",
                "class_a_ppo",
                "class_a_bismsb",
                "class_s_ppo",
                "class_s_bismsb"};
        }
    }
    else
    {
        if (!corr)
        {
            orderedNamesVec = new std::vector<std::string>{
                "deltam21",
                theta12name,
                "reactor_nubar",
                "geonu_U",
                "geonu_Th",
                "alphan_CScatter",
                "alphan_OExcited",
                "alphan_PRecoil",
                "bipolike",
                "atmospheric",
                "energy_scale",
                "energy_conv",
                "birks_constant",
                "p_recoil_energy_scale",
                "reactor_nubar2",
                "geonu_U2",
                "geonu_Th2",
                "alphan_CScatter2",
                "alphan_OExcited2",
                "alphan_PRecoil2",
                "bipolike2",
                "atmospheric2",
                "energy_scale2",
                "energy_conv2",
                "birks_constant2",
                "p_recoil_energy_scale2"};
        }
        else
        {
            orderedNamesVec = new std::vector<std::string>{
                "deltam21",
                theta12name,
                "reactor_nubar_norm",
                "geonu_U_norm",
                "geonu_Th_norm",
                "alphan_CScatter_norm",
                "alphan_OExcited_norm",
                "alphan_PRecoil_norm",
                "bipolike_norm",
                "atmospheric_norm",
                "energy_scale",
                "energy_conv",
                "birks_constant",
                "p_recoil_energy_scale",
                "energy_scale2",
                "energy_conv2",
                "birks_constant2",
                "p_recoil_energy_scale2"};
        }
    }

    std::vector<std::string> *selectedNames = new std::vector<std::string>;
    for (const std::string &name : *orderedNamesVec)
        if (std::find(namesVec->begin(), namesVec->end(), name) != namesVec->end())
            selectedNames->push_back(name);
    for (const std::string &name : *namesVec)
        if (std::find(selectedNames->begin(), selectedNames->end(), name) == selectedNames->end())
            selectedNames->push_back(name);
    orderedNamesVec = selectedNames;

    std::vector<double> *tempNomsVec = new std::vector<double>(nParams, 0.0);
    std::vector<double> *tempErrVec = new std::vector<double>(nParams, 0.0);
    std::vector<double> *tempConstrMeansVec = new std::vector<double>(nParams, 0.0);
    std::vector<double> *tempConstrErrsVec = new std::vector<double>(nParams, 0.0);
    std::vector<std::string> *tempLabelsVec = new std::vector<std::string>(nParams, "");
    TMatrixT<double> *tempCovMatrix = new TMatrixT<double>(nParams, nParams);

    // Map original parameter names to their old indices
    std::unordered_map<std::string, int> nameToOldIndex;
    for (int iParam = 0; iParam < nParams; ++iParam)
    {
        nameToOldIndex[namesVec->at(iParam)] = iParam;
    }

    // Now loop over the new vector (that's already in the correct order), and find the index of the same name
    // in the old vector, and set the noms and constraints to be the ones for that index
    for (int iParam = 0; iParam < nParams; ++iParam)
    {
        int iOriginal = nameToOldIndex[orderedNamesVec->at(iParam)];
        tempNomsVec->at(iParam) = nomsVec->at(iOriginal);
        tempErrVec->at(iParam) = errVec->at(iOriginal);
        tempConstrMeansVec->at(iParam) = constrMeansVec->at(iOriginal);
        tempConstrErrsVec->at(iParam) = constrErrsVec->at(iOriginal);
        tempLabelsVec->at(iParam) = labelsVec->at(iOriginal);

        for (int jParam = 0; jParam < nParams; ++jParam)
        {
            int jOriginal = nameToOldIndex[orderedNamesVec->at(jParam)];
            (*tempCovMatrix)(iParam, jParam) = (*covMatrix)(iOriginal, jOriginal);
        }
    }
    for (int iParam = 0; iParam < nParams; ++iParam)
    {
        for (int jParam = 0; jParam < nParams; ++jParam)
        {
            if (orderedNamesVec->at(jParam) == "geonu_U_norm" && orderedNamesVec->at(iParam) == "geonu_Th_norm")
                std::cout << "Geo correlation: " << (*tempCovMatrix)(iParam, jParam) / (sqrt((*tempCovMatrix)(iParam, iParam)) * sqrt((*tempCovMatrix)(jParam, jParam))) << std::endl;
            else if (orderedNamesVec->at(jParam) == "geonu_U" && orderedNamesVec->at(iParam) == "geonu_Th")
                std::cout << "PPO Geo correlation: " << (*tempCovMatrix)(iParam, jParam) / (sqrt((*tempCovMatrix)(iParam, iParam)) * sqrt((*tempCovMatrix)(jParam, jParam))) << std::endl;
            else if (orderedNamesVec->at(jParam) == "geonu_U2" && orderedNamesVec->at(iParam) == "geonu_Th2")
                std::cout << "bisMSB Geo correlation: " << (*tempCovMatrix)(iParam, jParam) / (sqrt((*tempCovMatrix)(iParam, iParam)) * sqrt((*tempCovMatrix)(jParam, jParam))) << std::endl;
        }
    }
    namesVec = orderedNamesVec;
    labelsVec = tempLabelsVec;
    nomsVec = tempNomsVec;
    errVec = tempErrVec;
    constrMeansVec = tempConstrMeansVec;
    constrErrsVec = tempConstrErrsVec;
    covMatrix = tempCovMatrix;
}

inline void DrawFitParams(const char *filename, std::map<std::string, double> &branchValues,
                         std::vector<std::string> *paramNames, std::vector<double> *paramErr,
                         std::vector<std::string> *labelsVec, std::vector<double> *nomVals,
                         std::vector<double> *constrMeans, std::vector<double> *constrErr,
                         TMatrixT<double> *fullCovMatrix, const std::string &theta12name,
                         bool corr, bool alphanclass,
                         const std::map<std::string, std::pair<double, double>> *intervalErrors = nullptr)
{
    // Reorder vectors to the order we want to plot them
    sortVectors(paramNames, paramErr, labelsVec, nomVals, constrMeans, constrErr, fullCovMatrix, theta12name, corr, alphanclass);

    TH1D *hNom = new TH1D("hNominal", "Relative Nominal Values", paramNames->size(), 0, paramNames->size() - 1);
    TH1D *hConstr = new TH1D("hConstr", "Constraints Relative to Nominals", paramNames->size(), 0, paramNames->size() - 1);
    TH1D *hPostFit = new TH1D("hPostFit", "Postfit Values Relative to Nominals", paramNames->size(), 0, paramNames->size() - 1);
    TGraphAsymmErrors *posterior = intervalErrors ? new TGraphAsymmErrors(paramNames->size()) : nullptr;

    for (int iParam = 0; iParam < paramNames->size(); iParam++)
    {
        // For reactor nubar, rate in fit config is unoscillated, but postfit value is oscillated
        // So we need to multiply by the ratio saved in outputted fit file
        if (paramNames->at(iParam) == reactorpar1)
        {
            nomVals->at(iParam) = nomVals->at(iParam) * branchValues[reactorpar1 + "_ratio"];
            constrMeans->at(iParam) = constrMeans->at(iParam) * branchValues[reactorpar1 + "_ratio"];
            constrErr->at(iParam) = constrErr->at(iParam) * branchValues[reactorpar1 + "_ratio"];
        }
        else if (paramNames->at(iParam) == reactorpar2)
        {
            nomVals->at(iParam) = nomVals->at(iParam) * branchValues[reactorpar2 + "_ratio"];
            constrMeans->at(iParam) = constrMeans->at(iParam) * branchValues[reactorpar2 + "_ratio"];
            constrErr->at(iParam) = constrErr->at(iParam) * branchValues[reactorpar2 + "_ratio"];
        }

        double nom = 1;
        double constrmean = 1;
        double postfit = 1;
        if (nomVals->at(iParam) == 0)
        {
            nom = 1.0;
            constrmean = constrMeans->at(iParam) + 1;
            postfit = branchValues[paramNames->at(iParam)] + 1;
        }
        else
        {
            nom = nomVals->at(iParam);
            constrmean = constrMeans->at(iParam);
            postfit = branchValues[paramNames->at(iParam)];
        }

        hNom->SetBinContent(iParam + 1, nom / nom);
        hConstr->GetXaxis()->SetBinLabel(iParam + 1, labelsVec->at(iParam).c_str());
        hConstr->SetBinContent(iParam + 1, constrmean / nom);
        hConstr->SetBinError(iParam + 1, constrErr->at(iParam) / nom);
        hPostFit->SetBinContent(iParam + 1, postfit / nom);
        hPostFit->SetBinError(iParam + 1, paramErr->at(iParam) / nom);
        if (posterior)
        {
            const std::pair<double, double> errors = intervalErrors->at(paramNames->at(iParam));
            posterior->SetPoint(iParam, hPostFit->GetBinCenter(iParam + 1), postfit / nom);
            posterior->SetPointError(iParam, 0, 0, errors.first / std::abs(nom), errors.second / std::abs(nom));
        }
        std::cout << "Par: " << paramNames->at(iParam) << std::endl;
        std::cout << "Nom Mean: " << nomVals->at(iParam) << std::endl;
        std::cout << "Constr: " << constrMeans->at(iParam) << " " << constrErr->at(iParam) << std::endl;
        std::cout << "Fit: " << branchValues[paramNames->at(iParam)] << " " << paramErr->at(iParam) << std::endl;
        if (paramNames->at(iParam) == reactorpar1)
            std::cout << "Osc Ratio: " << branchValues[reactorpar1 + "_ratio"] << std::endl;
        else if (paramNames->at(iParam) == reactorpar2)
        {
            std::cout << "Osc Ratio: " << branchValues[reactorpar2 + "_ratio"] << std::endl;
        }
        else if ((paramNames->at(iParam).find("reactor") != std::string::npos))
        {
            std::cout << "PPO Osc Ratio: " << branchValues[reactorpar1 + "_ratio"] << std::endl;
            std::cout << "bisMSB Osc Ratio: " << branchValues[reactorpar2 + "_ratio"] << std::endl;
        }
        std::cout << std::endl;
    }

    // Draw the histograms
    TCanvas *c1 = new TCanvas("c1", "Params", 1500, 1000);
    c1->SetBottomMargin(0.18);
    c1->SetRightMargin(0.12);
    gPad->SetFrameLineWidth(2);
    gStyle->SetOptStat(0);
    gPad->SetGrid(1);

    hNom->SetLineColor(kGreen);
    hNom->SetLineWidth(2);
    hConstr->SetLineColor(kBlue);
    hConstr->SetLineWidth(2);
    hPostFit->SetLineColor(kRed);
    hPostFit->SetLineWidth(2);

    hConstr->GetYaxis()->SetRangeUser(0, 2);
    if (posterior)
    {
        double lower = 0;
        double upper = 2;
        for (int i = 0; i < posterior->GetN(); ++i)
        {
            const double value = posterior->GetY()[i];
            lower = std::min(lower, value - posterior->GetErrorYlow(i));
            upper = std::max(upper, value + posterior->GetErrorYhigh(i));
        }
        hConstr->GetYaxis()->SetRangeUser(lower, upper + 0.05 * (upper - lower));
    }
    hConstr->GetYaxis()->SetTitle("Relative to Nominal");
    hConstr->GetYaxis()->SetTitleOffset(1.2);
    hConstr->GetXaxis()->SetLabelOffset(0.007);
    hConstr->GetXaxis()->SetTitle("Fit Parameters");
    hConstr->GetXaxis()->SetTitleOffset(2.2);
    hConstr->SetTitle("");

    hConstr->GetXaxis()->SetTitleFont(42);
    hConstr->GetYaxis()->SetTitleFont(42);
    hConstr->GetXaxis()->SetLabelFont(42);
    hConstr->GetYaxis()->SetLabelFont(42);
    hConstr->GetXaxis()->SetLabelSize(0.03);
    hConstr->SetTitleFont(42);

    hConstr->Draw("E1");
    if (posterior)
    {
        posterior->SetLineColor(kRed);
        posterior->SetLineWidth(2);
        posterior->SetMarkerStyle(20);
        posterior->SetMarkerColor(kRed);
        posterior->Draw("PZ SAME");
    }
    else
        hPostFit->Draw("E1same");

    TLegend *t1 = new TLegend(0.65, 0.73, 0.85, 0.88);
    t1->AddEntry(hConstr, "Prefit", "l");
    if (posterior)
        t1->AddEntry(posterior, "MCMC median, 68% interval", "lp");
    else
        t1->AddEntry(hPostFit, "Postfit", "l");
    t1->SetLineWidth(2);
    t1->SetTextFont(42);
    t1->Draw();

    // Save plot as image and rootfile
    struct stat st = {0};
    std::filesystem::path pathObj(filename);
    pathObj.replace_filename("plots/");
    if (stat(pathObj.string().c_str(), &st) == -1)
        mkdir(pathObj.string().c_str(), 0700);
    pathObj.replace_filename("params.pdf");
    c1->SaveAs(pathObj.string().c_str());
    pathObj.replace_filename("params.root");
    TFile *outfile = new TFile(pathObj.string().c_str(), "RECREATE");
    outfile->cd();
    hNom->Write("nominal");
    hConstr->Write("constraints");
    hPostFit->Write("postfit");
    if (posterior)
    {
        posterior->Write("posterior_intervals");
        TNamed("uncertainty_method", "posterior_intervals: median and equal-tailed 68% interval; postfit histogram errors: posterior standard deviation").Write();
    }
    c1->Write("c1");

    // Now plot the correlation matrix
    int nParams = fullCovMatrix->GetNrows();
    TH2D *hCorrMatrix = new TH2D("hCorrMatrix", "Correlation Matrix", nParams, 0, nParams, nParams, 0, nParams);

    for (int i = 0; i < nParams; ++i)
    {
        for (int j = 0; j < nParams; ++j)
        {
            hCorrMatrix->SetBinContent(i + 1, j + 1, (*fullCovMatrix)(i, j) / (sqrt((*fullCovMatrix)(i, i)) * sqrt((*fullCovMatrix)(j, j))));
        }
    }

    // Set axis labels
    for (int i = 0; i < nParams; ++i)
    {
        hCorrMatrix->GetXaxis()->SetBinLabel(i + 1, labelsVec->at(i).c_str());
        hCorrMatrix->GetYaxis()->SetBinLabel(i + 1, labelsVec->at(i).c_str());
    }

    // Define a nice red to blue palette
    const Int_t NRGBs = 5;
    const Int_t NCont = 255;

    // White in the middle, darker red and blue at ends
    Double_t stops[NRGBs] = {0.00, 0.25, 0.50, 0.75, 1.00};
    Double_t red[NRGBs] = {0.00, 0.25, 1.00, 1.00, 0.50};
    Double_t green[NRGBs] = {0.00, 0.25, 1.00, 0.25, 0.00};
    Double_t blue[NRGBs] = {0.50, 1.00, 1.00, 0.25, 0.00};

    TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
    gStyle->SetNumberContours(NCont);

    // And make the plot
    TCanvas *c2 = new TCanvas("c2", "Correlations", 1500, 1000);
    c2->SetBottomMargin(0.18);
    c2->SetRightMargin(0.16);
    c2->SetLeftMargin(0.16);
    gPad->SetFrameLineWidth(2);
    gStyle->SetOptStat(0);
    gPad->SetGrid(1);

    hCorrMatrix->GetXaxis()->SetLabelOffset(0.007);
    hCorrMatrix->GetXaxis()->SetLabelFont(42);
    hCorrMatrix->GetYaxis()->SetLabelFont(42);
    hCorrMatrix->GetZaxis()->SetLabelFont(42);
    hCorrMatrix->GetZaxis()->SetRangeUser(-1, 1);
    hCorrMatrix->GetXaxis()->SetLabelSize(0.03);
    hCorrMatrix->GetYaxis()->SetLabelSize(0.03);
    hCorrMatrix->SetTitle("");
    hCorrMatrix->Draw("colz");

    outfile->cd();
    hCorrMatrix->Write("correlations");
    c2->Write("c2");
    pathObj.replace_filename("correlations.pdf");
    c2->SaveAs(pathObj.string().c_str());
}

}
