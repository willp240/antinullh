// ROOT Headers
#include <TFile.h>
#include <TTree.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TMatrixT.h>
#include <TParameter.h>
#include "PlotFitParams.hh"

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

void plotFixedOscParams(const char *filename = "fit_results.root", bool corr = true, bool alphanclass = false)
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

    // Floated fits already contain oscillation errors and the full Minuit covariance
    TParameter<int> *floatedFlag = nullptr;
    file->GetObject("oscillation_floated", floatedFlag);
    TMatrixT<double> *fullCovMatrix = covMatrix;
    if (!floatedFlag || !floatedFlag->GetVal())
    {
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
        fullCovMatrix = appendVectorsToMatrix(covMatrix, dmCovVec, thCovVec);

    }

    antinupostfit::DrawFitParams(filename, branchValues, paramNames, paramErr, labelsVec,
                               nomVals, constrMeans, constrErr, fullCovMatrix, theta12name,
                               corr, alphanclass);
}
