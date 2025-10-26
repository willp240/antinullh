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

    int nParams = namesVec->size();
    // This is the order we want to plot them in (osc, signal, geo, alpha n, other bgs, systematics)
    /*std::vector<std::string> *orderedNamesVec = new std::vector<std::string>{
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
        "class_s_bismsb"
    };*/

    std::vector<std::string> *orderedNamesVec = new std::vector<std::string>{
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
        "class_s_bismsb"
    };

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

    namesVec = orderedNamesVec;
    labelsVec = tempLabelsVec;
    nomsVec = tempNomsVec;
    errVec = tempErrVec;
    constrMeansVec = tempConstrMeansVec;
    constrErrsVec = tempConstrErrsVec;
    covMatrix = tempCovMatrix;
}

TMatrixT<double> *appendVectorsToMatrix(const TMatrixT<double> *covMatrix, const std::vector<double> &v1, const std::vector<double> &v2)
{
    const int N = covMatrix->GetNrows();
    if (covMatrix->GetNcols() != N)
        throw std::runtime_error("Matrix must be square");
    if (v1.size() != static_cast<size_t>(N + 2) || v2.size() != static_cast<size_t>(N + 2))
        throw std::runtime_error("Vectors to be appended to NxN matrix must have size N+2");

    TMatrixT<double> *newCovMatrix = new TMatrixT<double>(N + 2, N + 2);
    // Top-left block: copy of the original matrix
    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j)
            (*newCovMatrix)(i, j) = (*covMatrix)(i, j);

    // New columns (N and N+1) for rows 0..N-1
    for (int i = 0; i < N; ++i)
    {
        (*newCovMatrix)(i, N) = v1[i];
        (*newCovMatrix)(i, N + 1) = v2[i];
    }

    // New rows (N and N+1) for cols 0..N-1 (mirror to keep symmetry)
    for (int j = 0; j < N; ++j)
    {
        (*newCovMatrix)(N, j) = v1[j];
        (*newCovMatrix)(N + 1, j) = v2[j];
    }

    // Bottom-right 2x2 from the tails of v1, v2
    (*newCovMatrix)(N, N) = v1[N];
    (*newCovMatrix)(N, N + 1) = v1[N + 1];
    (*newCovMatrix)(N + 1, N) = v2[N];
    (*newCovMatrix)(N + 1, N + 1) = v2[N + 1];

    return newCovMatrix;
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
    std::vector<double> *paramVals = nullptr;
    std::vector<double> *paramErr = nullptr;
    std::vector<double> *constrMeans = nullptr;
    std::vector<double> *constrErr = nullptr;
    std::vector<std::string> *labelsVec = nullptr;
    TMatrixT<double> *covMatrix = nullptr;

    file->GetObject("param_names", paramNames);
    file->GetObject("param_nom_values", nomVals);
    file->GetObject("param_fit_values", paramVals);
    file->GetObject("param_fit_err", paramErr);
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

    // Now put the oscillation errors into the error vector. This will still be in the same order as the param values vector etc.
    paramErr->push_back(std::max(dmLowErr, dmHighErr));
    paramErr->push_back(std::max(thLowErr, thHighErr));

    // Now we'll find covariance of each parameter with each osc par and put in vectors
    std::vector<double> dmCovVec;
    std::vector<double> thCovVec;

    double mindeltam = tree->GetMinimum("deltam21");
    double maxdeltam = tree->GetMaximum("deltam21");
    double mintheta = tree->GetMinimum(theta12name.c_str());
    double maxtheta = tree->GetMaximum(theta12name.c_str());

    for (int iPar = 0; iPar < paramNames->size(); iPar++)
    {

        double parmin = tree->GetMinimum(paramNames->at(iPar).c_str());
        double parmax = tree->GetMaximum(paramNames->at(iPar).c_str());

        TH2D hdeltam("hdeltamcov", ("deltam21;" + paramNames->at(iPar)).c_str(), 500, mindeltam, maxdeltam, 500, parmin, parmax);
        std::string dmDrawCmd = paramNames->at(iPar) + ":deltam21" + " >> hdeltamcov";
        tree->Draw(dmDrawCmd.c_str(), "LLH", "goff");
        double deltamCov = hdeltam.GetCovariance();
        dmCovVec.push_back(deltamCov);

        TH2D htheta("hthetacov", (theta12name + ";" + paramNames->at(iPar)).c_str(), 500, mintheta, maxtheta, 500, parmin, parmax);
        std::string thDrawCmd = paramNames->at(iPar) + ":" + theta12name + " >> hthetacov";
        tree->Draw(thDrawCmd.c_str(), "LLH", "goff");
        double thetaCov = htheta.GetCovariance();
        thCovVec.push_back(thetaCov);
    }

    tree->GetEntry(bestLLHEntry);

    // Now append these oscillation covariance vectors to the covariance matrix
    TMatrixT<double> *fullCovMatrix = appendVectorsToMatrix(covMatrix, dmCovVec, thCovVec);

    // Reorder vectors to the order we want to plot them
    sortVectors(paramNames, paramErr, labelsVec, nomVals, constrMeans, constrErr, fullCovMatrix, theta12name);

    TH1D *hNom = new TH1D("hNominal", "Relative Nominal Values", paramNames->size(), 0, paramNames->size() - 1);
    TH1D *hConstr = new TH1D("hConstr", "Constraints Relative to Nominals", paramNames->size(), 0, paramNames->size() - 1);
    TH1D *hPostFit = new TH1D("hPostFit", "Postfit Values Relative to Nominals", paramNames->size(), 0, paramNames->size() - 1);

    for (int iParam = 0; iParam < paramNames->size(); iParam++)
    {
        // For reactor nubar, rate in fit config is unoscillated, but postfit value is oscillated.
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
        if(nomVals->at(iParam) == 0)
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
        std::cout << "Par: " << paramNames->at(iParam) << std::endl;
        std::cout << "Nom Mean: " << nomVals->at(iParam) << std::endl;
        std::cout << "Constr: " << constrMeans->at(iParam) << " " << constrErr->at(iParam) << std::endl;
        std::cout << "Fit: " << branchValues[paramNames->at(iParam)] << " " << paramErr->at(iParam) << std::endl
                  << std::endl;
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
    c2->SetRightMargin(0.155);
    c2->SetLeftMargin(0.15);
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
