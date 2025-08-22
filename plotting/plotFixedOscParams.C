// ROOT Headers
#include <TFile.h>
#include <TTree.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TMatrixT.h>

// c++ Headers
#include <sys/stat.h>
#include <unordered_map>

/* ///////////////////////////////////////////////////////////////////
///
/// Script for plotting post fit parameter values and prefit
/// constraints relative to nominal values for fixed oscillation fits,
/// as well as the postfit covariance matrix.
///
/// The user inputs the root file made by makeFixedOscTree, and it
/// first finds the entry with the lowest LLH. All the branches for
/// this entry are read into a map, and the nominal values and
/// constraints vectors are read in from the input file too. Those
/// vectors are sorted into the order we want to plot (so this will
/// need to be updated when parameters change). The covariance matrix
/// is loaded and the entries rearranged to the order we want to plot
/// them in, and the correlations are calculated and plotted.
///
/// The plots are drawn and the canvases are saved as pdf and root files,
/// along with each of the histograms.
///
/////////////////////////////////////////////////////////////////// */

// Hard code reactor parameter names to make sure we use the right ratio for each
// Inelegant, but we hard code the parameter order anyway
const std::string reactorpar1 = "reactor_nubar";
const std::string reactorpar2 = "reactor_nubar2";

/// Function to sort vectors of names, nominals, and constraints into the order we want to plot them in.
void sortVectors(std::vector<std::string> *&namesVec, std::vector<double> *&errVec, std::vector<std::string> *&labelsVec, std::vector<double> *&nomsVec,
                 std::vector<double> *&constrMeansVec, std::vector<double> *&constrErrsVec, TMatrixT<double> *&covMatrix, std::string theta12name)
{

    // This is the order we want to plot them in (osc, signal, geo, alpha n, other bgs, systematics)
    std::vector<std::string> *tempNamesVec = new std::vector<std::string>{"deltam21",
                                                                          theta12name,
                                                                          "reactor_nubar",
                                                                          "reactor_nubar2",
                                                                          "geonu_U",
                                                                          "geonu_U2",
                                                                          "geonu_Th",
                                                                          "geonu_Th2",
                                                                          "alphan_CScatter",
                                                                          "alphan_CScatter2",
                                                                          "alphan_OExcited",
                                                                          "alphan_OExcited2",
                                                                          "alphan_PRecoil",
                                                                          "alphan_PRecoil2",
                                                                          "energy_scale",
                                                                          "energy_scale2",
                                                                          "energy_conv",
                                                                          "energy_conv2",
                                                                          "birks_constant",
                                                                          "birks_constant2",
                                                                          "p_recoil_energy_scale",
                                                                          "p_recoil_energy_scale2"};

    int nParams = namesVec->size();

    std::vector<double> *tempNomsVec = new std::vector<double>(nParams, 0.0);
    std::vector<double> *tempErrVec = new std::vector<double>(nParams, 0.0);
    std::vector<double> *tempConstrMeansVec = new std::vector<double>(nParams, 0.0);
    std::vector<double> *tempConstrErrsVec = new std::vector<double>(nParams, 0.0);
    std::vector<std::string> *tempLabelsVec = new std::vector<std::string>(nParams, "");
    // For the fixed osc fits, the matrix won't have theta12 and deltam
    TMatrixT<double> *tempCovMatrix = new TMatrixT<double>(nParams - 2, nParams - 2);
    std::vector<std::string> *matrixNamesVec = new std::vector<std::string>(namesVec->begin(), namesVec->end());

    // Matrix names won't include osc params
    std::vector<std::string> *tempMatrixNamesVec = new std::vector<std::string>(tempNamesVec->begin() + 2, tempNamesVec->end());

    // Map original parameter names to their new indices
    std::unordered_map<std::string, int> nameToOldIndex;
    for (int iParam = 0; iParam < nParams; ++iParam)
    {
        nameToOldIndex[namesVec->at(iParam)] = iParam;
    }

    // Now loop over the new vector (that's already in the correct order), and find the index of the same name
    // in the old vector, and set the noms and constraints to be the ones for that index
    for (int iParam = 0; iParam < nParams; ++iParam)
    {
        int iOriginal = nameToOldIndex[tempNamesVec->at(iParam)];

        tempNomsVec->at(iParam) = nomsVec->at(iOriginal);
        tempErrVec->at(iParam) = errVec->at(iOriginal);
        tempConstrMeansVec->at(iParam) = constrMeansVec->at(iOriginal);
        tempConstrErrsVec->at(iParam) = constrErrsVec->at(iOriginal);
        tempLabelsVec->at(iParam) = labelsVec->at(iOriginal);
    }

    // And redo this for the matrix indices which will be different (no theta or deltam)
    for (int iParam = 0; iParam < matrixNamesVec->size(); ++iParam)
    {
        // Erase the deltam and theta entries from the matrix name vector to get the right order for those later
        if (matrixNamesVec->at(iParam) == "deltam21" || matrixNamesVec->at(iParam) == theta12name)
        {
            matrixNamesVec->erase(matrixNamesVec->begin() + iParam) - 1;
            iParam--;
        }
    }

    int nMatrixParams = matrixNamesVec->size();

    // Map original parameter names to their new indices
    std::unordered_map<std::string, int> nameToOldMatrixIndex;
    for (int iParam = 0; iParam < nMatrixParams; ++iParam)
    {
        nameToOldMatrixIndex[matrixNamesVec->at(iParam)] = iParam;
    }

    for (int iParam = 0; iParam < nMatrixParams; ++iParam)
    {
        for (int jParam = 0; jParam < nMatrixParams; ++jParam)
        {
            int iOriginal = nameToOldMatrixIndex[tempMatrixNamesVec->at(iParam)];
            int jOriginal = nameToOldMatrixIndex[tempMatrixNamesVec->at(jParam)];
            (*tempCovMatrix)(iParam, jParam) = (*covMatrix)(iOriginal, jOriginal);
        }
    }

    namesVec = tempNamesVec;
    labelsVec = tempLabelsVec;
    nomsVec = tempNomsVec;
    errVec = tempErrVec;
    constrMeansVec = tempConstrMeansVec;
    constrErrsVec = tempConstrErrsVec;
    covMatrix = tempCovMatrix;
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

    // Check the form of theta 12
    std::string theta12name;
    if (tree->GetBranch("theta12"))
        theta12name = "theta12";
    else if (tree->GetBranch("sintheta12"))
        theta12name = "sintheta12";
    else if (tree->GetBranch("sinsqtheta12"))
        theta12name = "sinsqtheta12";

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
    std::cout << "Best theta: " << branchValues[theta12name] << std::endl;
    std::cout << "Best deltam: " << branchValues["deltam21"] << std::endl
              << std::endl;

    // Retrieve the vectors from the file: nominal, constr means, constr err, constr names, labels, cov matrix
    std::vector<std::string> *paramNames = nullptr;
    std::vector<double> *nomVals = nullptr;
    std::vector<double> *constrMeans = nullptr;
    std::vector<double> *constrErr = nullptr;
    std::vector<std::string> *labelsVec = nullptr;
    std::vector<double> *paramErr = nullptr;
    TMatrixT<double> *covMatrix = nullptr;

    file->GetObject("param_names", paramNames);
    file->GetObject("param_asimov_values", nomVals);
    file->GetObject("constr_mean_values", constrMeans);
    file->GetObject("constr_sigma_values", constrErr);
    file->GetObject("tex_labels", labelsVec);
    file->GetObject("covMatrix", covMatrix);

    // Open the llh plots file if it exists to get the oscillation uncertainties
    std::filesystem::path filePath(filename);
    std::filesystem::path fileDir = filePath.parent_path();
    std::ostringstream llhFileName;

    llhFileName << fileDir.string() << "/plots/LLH.root";

    TFile *llhFile = new TFile(llhFileName.str().c_str(), "READ");
    double thLowErr = -999;
    double thHighErr = -999;
    double dmLowErr = -999;
    double dmHighErr = -999;

    if (!llhFile || llhFile->IsZombie())
    {
        std::cerr << "Error: Could not open file " << llhFileName.str()
                  << ". Defaulting to grid space oscillation parameter errors" << std::endl;
        thLowErr = 0.18;
        thHighErr = 0.18;
        dmLowErr = 7E-8;
        dmHighErr = 7E-8;
    }
    else
    {

        std::map<std::string, double> *dmSigmas = llhFile->Get<std::map<std::string, double>>("dmSigmas");
        std::map<std::string, double> *thSigmas = llhFile->Get<std::map<std::string, double>>("thSigmas");

        if (!dmSigmas || !thSigmas)
        {
            std::cerr << "Error: Could not find maps in file " << llhFileName.str()
                      << ". Defaulting to grid space oscillation parameter errors" << std::endl;
            if (theta12name == "theta12")
            {
                thLowErr = 0.18;
                thHighErr = 0.18;
            }
            else if (theta12name == "sintheta12" || theta12name == "sinsqtheta12")
            {
                thLowErr = 0.002;
                thHighErr = 0.1;
            }
            dmLowErr = 7E-8;
            dmHighErr = 7E-8;
        }
        else
        {
            thLowErr = abs((*thSigmas)["left"] - (*thSigmas)["bestfit"]);
            thHighErr = abs((*thSigmas)["right"] - (*thSigmas)["bestfit"]);
            dmLowErr = abs((*dmSigmas)["left"] - (*dmSigmas)["bestfit"]);
            dmHighErr = abs((*dmSigmas)["right"] - (*dmSigmas)["bestfit"]);
        }
    }

    std::ostringstream fitFileName;

    fitFileName << fileDir.string() << "/th" << std::fixed << std::setprecision(3) << branchValues[theta12name] << "/th"
                << std::fixed << std::setprecision(3) << branchValues[theta12name] << "_dm"
                << std::fixed << std::setprecision(8) << branchValues["deltam21"] << "/fit_result.root";

    TFile *fitFile = new TFile(fitFileName.str().c_str(), "READ");
    fitFile->GetObject("paramErr", paramErr);

    // Now find the position of deltam and theta in the other vectors, so we can insert them into paramErr at the same place
    // before reordering
    std::vector<std::string>::iterator deltampos = std::find(paramNames->begin(), paramNames->end(), "deltam21");
    std::vector<std::string>::iterator thetapos = std::find(paramNames->begin(), paramNames->end(), theta12name);
    paramErr->insert(paramErr->begin() + std::distance(paramNames->begin(), deltampos), std::max(dmLowErr, dmHighErr));
    paramErr->insert(paramErr->begin() + std::distance(paramNames->begin(), thetapos), std::max(thLowErr, thHighErr));

    // Reorder vectors to the order we want to plot them
    sortVectors(paramNames, paramErr, labelsVec, nomVals, constrMeans, constrErr, covMatrix, theta12name);

    TH1D *hNom = new TH1D("hNominal", "Relative Nominal Values", paramNames->size(), 0, paramNames->size() - 1);
    TH1D *hConstr = new TH1D("hConstr", "Constraints Relative to Nominals", paramNames->size(), 0, paramNames->size() - 1);
    TH1D *hPostFit = new TH1D("hPostFit", "Postfit Values Relative to Nominals", paramNames->size(), 0, paramNames->size() - 1);

    for (int iParam = 0; iParam < paramNames->size(); iParam++)
    {
        // For reactor nubar, rate in fit config is unoscillated, but postfit value is oscillated.
        // So we need to multiply by the ratio saved in outputted fit file
        if (paramNames->at(iParam).find(reactorpar1) != std::string::npos)
        {
            nomVals->at(iParam) = nomVals->at(iParam) * branchValues[reactorpar1 + "_ratio"];
            constrMeans->at(iParam) = constrMeans->at(iParam) * branchValues[reactorpar1 + "_ratio"];
            constrErr->at(iParam) = constrErr->at(iParam) * branchValues[reactorpar1 + "_ratio"];
        }
        else if (paramNames->at(iParam).find(reactorpar2) != std::string::npos)
        {
            nomVals->at(iParam) = nomVals->at(iParam) * branchValues[reactorpar2 + "_ratio"];
            constrMeans->at(iParam) = constrMeans->at(iParam) * branchValues[reactorpar2 + "_ratio"];
            constrErr->at(iParam) = constrErr->at(iParam) * branchValues[reactorpar2 + "_ratio"];
        }

        hNom->SetBinContent(iParam + 1, nomVals->at(iParam) / nomVals->at(iParam));
        hConstr->GetXaxis()->SetBinLabel(iParam + 1, labelsVec->at(iParam).c_str());
        hConstr->SetBinContent(iParam + 1, constrMeans->at(iParam) / nomVals->at(iParam));
        hConstr->SetBinError(iParam + 1, constrErr->at(iParam) / nomVals->at(iParam));
        hPostFit->SetBinContent(iParam + 1, branchValues[paramNames->at(iParam)] / nomVals->at(iParam));
        hPostFit->SetBinError(iParam + 1, paramErr->at(iParam) / nomVals->at(iParam));
        std::cout << "Par: " << paramNames->at(iParam) << std::endl;
        std::cout << "Nom Mean: " << nomVals->at(iParam) << std::endl;
        std::cout << "Constr: " << constrMeans->at(iParam) << " " << constrErr->at(iParam) << std::endl;
        std::cout << "Fit: " << branchValues[paramNames->at(iParam)] << " " << paramErr->at(iParam) << std::endl
                  << std::endl;
    }

    // Draw the histograms
    TCanvas *c1 = new TCanvas("c1", "Params", 800, 600);
    c1->SetBottomMargin(0.18);
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
    hConstr->GetYaxis()->SetTitle("Relative to Nominal");
    hConstr->GetYaxis()->SetTitleOffset(1.2);
    hConstr->GetXaxis()->SetLabelOffset(0.007);
    hConstr->GetXaxis()->SetTitle("Fit Parameters");
    hConstr->GetXaxis()->SetTitleOffset(2.0);
    hConstr->SetTitle("");

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
    c1->Write("c1");

    // Now plot the correlation matrix
    int nParams = covMatrix->GetNrows();
    TH2D *hCorrMatrix = new TH2D("hCorrMatrix", "Correlation Matrix", nParams, 0, nParams, nParams, 0, nParams);

    for (int i = 0; i < nParams; ++i)
    {
        for (int j = 0; j < nParams; ++j)
        {
            hCorrMatrix->SetBinContent(i + 1, j + 1, (*covMatrix)(i, j) / (sqrt((*covMatrix)(i, i)) * sqrt((*covMatrix)(j, j))));
        }
    }

    // Set axis labels
    for (int i = 0; i < nParams; ++i)
    {
        hCorrMatrix->GetXaxis()->SetBinLabel(i + 1, labelsVec->at(i + 2).c_str());
        hCorrMatrix->GetYaxis()->SetBinLabel(i + 1, labelsVec->at(i + 2).c_str());
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
    TCanvas *c2 = new TCanvas("c2", "Correlations", 800, 600);
    c1->SetBottomMargin(0.18);
    c2->SetRightMargin(0.15);
    c2->SetLeftMargin(0.15);
    gPad->SetFrameLineWidth(2);
    gStyle->SetOptStat(0);
    gPad->SetGrid(1);

    hCorrMatrix->GetXaxis()->SetLabelOffset(0.007);
    hCorrMatrix->GetXaxis()->SetLabelFont(42);
    hCorrMatrix->GetYaxis()->SetLabelFont(42);
    hCorrMatrix->GetZaxis()->SetLabelFont(42);
    hCorrMatrix->GetZaxis()->SetRangeUser(-1, 1);
    hCorrMatrix->SetTitle("");
    hCorrMatrix->Draw("colz");

    outfile->cd();
    hCorrMatrix->Write("correlations");
    c2->Write("c2");
    pathObj.replace_filename("correlations.pdf");
    c2->SaveAs(pathObj.string().c_str());
}
