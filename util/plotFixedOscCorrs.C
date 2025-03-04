#include <TFile.h>
#include <TTree.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TStyle.h>

/* ///////////////////////////////////////////////////////////////////
///

///
/////////////////////////////////////////////////////////////////// */

/// Function to sort vectors of names and labels into the order we want them in the matrix plots
void makeIndexMap(std::vector<std::string> *namesVec, std::vector<std::string> *labelsVec, std::map<int, int> &indexMap)
{

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
    // in the old vector
    for (int iTempName = 2; iTempName < tempNamesVec->size(); iTempName++)
    {
        int offset = 0;
        for (int iName = 0; iName < namesVec->size(); iName++)
        {
            if (namesVec->at(iName) == tempNamesVec->at(0) || namesVec->at(iName) == tempNamesVec->at(1))
                offset++;
            if (tempNamesVec->at(iTempName) == namesVec->at(iName))
            {
                indexMap[iName - offset] = iTempName;
                std::cout << iName - offset << " " << iTempName << " " << offset << std::endl;
                break;
            }
        }
    }

    // Now check there aren't any parameters we weren't expecting. If there are, add them to the end
    int offset = 0;
    for (int iName = 0; iName < namesVec->size(); iName++)
    {
        if (namesVec->at(iName) == tempNamesVec->at(0) || namesVec->at(iName) == tempNamesVec->at(1))
            offset++;
        if (!indexMap[iName - offset] && namesVec->at(iName) != tempNamesVec->at(0) && namesVec->at(iName) != tempNamesVec->at(1))
        {
            indexMap[iName - offset] = indexMap.size();
        }
    }
}

void plotFixedOscCorrs(const char *filename = "fit_results.root")
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

    // Retrieve the vectors from the file: names, labels
    std::vector<std::string> *paramNames = nullptr;
    std::vector<std::string> *labelsVec = nullptr;

    file->GetObject("param_names", paramNames);
    file->GetObject("tex_labels", labelsVec);

    std::filesystem::path filePath(filename);
    std::filesystem::path fileDir = filePath.parent_path();
    std::ostringstream fitFileName;

    fitFileName << fileDir.string() << "/th" << std::fixed << std::setprecision(2) << branchValues["theta12"] << "/th"
                << std::fixed << std::setprecision(2) << branchValues["theta12"] << "_dm"
                << std::fixed << std::setprecision(8) << branchValues["deltam21"] << "/fit_result.root";

    TFile *fitFile = new TFile(fitFileName.str().c_str(), "OPEN");
    if (!fitFile->IsOpen())
    {
        std::cerr << "Error: Could not open file " << fitFileName.str() << std::endl;
        return;
    }

    // Load the TMatrixT<double>
    TMatrixT<double> *covMatrix = (TMatrixT<double> *)fitFile->Get("covMatrix");
    if (!covMatrix)
    {
        std::cerr << "Error: Could not find covMatrix in file " << fitFileName.str() << std::endl;
        return;
    }

    int N = covMatrix->GetNrows(); // Plus two for the oscillation parameters

    TH2D *hist = new TH2D("correlation_matrix", "Reordered Correlation Matrix",
                          N + 2, 0, N + 2, N + 2, 0, N + 2);


    // From here it should be:
    
    // Map labels/names vec elements to hist bins

    // Map covariance element to hist bin (using names vec as intermediary)

    // Loop over matrix elements and fill hist

    // Loop bins and set labels

    // Palette

    // Label and axis fonts, margins, titles

    // Save in root and pdf



    std::map<int, int> indexMap;
    makeIndexMap(paramNames, labelsVec, indexMap);

    // Loop over the matrix and reorder using the index map
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            int new_i = indexMap[i];
            int new_j = indexMap[j];
            // if(new_i < 2 || new_j < 2)
            //   continue;
            // std::cout << i << " " << new_i << std::endl;
            double corr = (*covMatrix)(i, j) / sqrt((*covMatrix)(i, i) * (*covMatrix)(j, j));
            std::cout << i << " " << " " << j << " " << (*covMatrix)(i, j) << std::endl;
            hist->SetBinContent(new_i + 1, new_j + 1, corr);
        }
    }

    for (int i = 0; i < hist->GetXaxis()->GetNbins(); i++)
    {
        hist->GetXaxis()->SetBinLabel(i+1, labelsVec->at(indexMap[i]).c_str());
        hist->GetYaxis()->SetBinLabel(i+1, labelsVec->at(indexMap[i]).c_str());
    }

    // Set color palette (red -> white -> blue)
    const Int_t NRGBs = 3;
    const Int_t NCont = 255;
    Double_t stops[NRGBs] = {0.0, 0.5, 1.0}; // Positions: red at -1, white at 0, blue at 1
    Double_t red[NRGBs] = {1.0, 1.0, 0.0};
    Double_t green[NRGBs] = {0.0, 1.0, 0.0};
    Double_t blue[NRGBs] = {0.0, 1.0, 1.0};
    TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
    gStyle->SetNumberContours(NCont);

    // Draw histogram
    TCanvas *c1 = new TCanvas("c1", "Reordered Correlation Matrix", 800, 600);
    hist->SetStats(0);
    hist->SetMinimum(-1); // Set min/max for color range
    hist->SetMaximum(1);
    hist->Draw("COLZ");

    // Save plot
    c1->SaveAs("correlation_matrix.png");

    // Clean up
    // fitFile->Close();

    /*
        // Save plot as image and rootfile
    std::filesystem::path pathObj(filename);
    pathObj.replace_filename("params.pdf");
    c1->SaveAs(pathObj.string().c_str());
    pathObj.replace_filename("params.root");
    TFile *outfile = new TFile(pathObj.string().c_str(), "RECREATE");
    outfile->cd();
    hNom->Write("nominal");
    hConstr->Write("constraints");
    hPostFit->Write("postfit");
    c1->Write("c1");*/
}
