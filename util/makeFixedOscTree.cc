// Antinu headers
#include <PDFConfigLoader.hh>
#include <DistBuilder.hh>
#include <FitConfigLoader.hh>
#include <EventConfigLoader.hh>
#include <SystConfigLoader.hh>
#include <SystFactory.hh>
#include <OscGridConfigLoader.hh>

// OXO headers
#include <ROOTNtuple.h>
#include <BinnedNLLH.h>
#include <IO.h>
#include <Rand.h>
#include <Minuit.h>

// ROOT headers
#include <TH1D.h>
#include <TFile.h>
#include <TTree.h>

// c++ headers
#include <sys/stat.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <string>

using namespace antinufit;

/* ///////////////////////////////////////////////////////////////////
///
/// App for making a TTree of the results of a set of fixed 
/// fits. Each Entry in the TTree represents one of the fits, and
/// each branch represents one of the fit parameters.
/// First open the root file (made by makeFixedOscTree) and find the
/// minimum LLH entry. The user supplies the fit config and 
/// oscgrid config used to produce one (any) of the fits.
///
/// There is a function for reading either the txt or root files
/// outputted by the fits. The root file is used by default as it
/// contains more information.
///
/// As well as the TTree with postfit values, the prefit values and
/// constraints are also saved in a file saved in the top directory
/// of the fit results
///
/////////////////////////////////////////////////////////////////// */

// Function to turn a map of string to double into two vectors, one of strings, one of doubles
// This is useful for writing to root files
void mapToVectors(std::map<std::string, double> sdmap, std::vector<std::string> &namesvec, std::vector<double> &meansvec)
{

    for (std::map<std::string, double>::iterator it = sdmap.begin(); it != sdmap.end(); ++it)
    {
        namesvec.push_back(it->first);
        meansvec.push_back(it->second);
    }
}

// Function to parse an outputted txt file and extract the fit parameters to a map
bool parseFitResultTxt(const std::string &filename, std::map<std::string, double> &branchMap)
{

    std::ifstream infile(filename);
    if (!infile)
    {
        std::cerr << "Failed to open file: " << filename << std::endl;
        return false;
    }

    std::string line;
    while (std::getline(infile, line))
    {
        std::istringstream iss(line);
        std::string key;
        double value;

        if (iss >> key >> value)
        {
            branchMap[key] = value;
            if (key == "LLH:")
                branchMap["LLH"] = value;
            if (key == "Fit Valid:")
                branchMap["fit_valid"] = value;
        }
    }

    return true;
}

// Function to parse an outputted root file and extract the fit parameters to a map
bool parseFitResultRoot(const std::string &filename, std::map<std::string, double> &branchMap)
{

    // Open the ROOT file
    TFile file(filename.c_str(), "READ");
    if (file.IsZombie())
    {
        std::cerr << "Error: Cannot open file " << filename << std::endl;
        return false;
    }

    // Retrieve the vectors from the file
    std::vector<std::string> *paramNames = nullptr;
    std::vector<double> *paramVals = nullptr;
    std::vector<double> *paramErr = nullptr;
    std::vector<double> *reacRatio = nullptr;

    file.GetObject("paramNames", paramNames);
    file.GetObject("paramVals", paramVals);
    file.GetObject("paramErr", paramErr);
    file.GetObject("reactorRatio", reacRatio);

    // Check if vectors were loaded correctly
    if (!paramNames || !paramVals || !paramErr)
    {
        std::cerr << "Error: One or more vectors could not be found in the file." << std::endl;
        return false;
    }

    if (paramNames->size() != paramVals->size())
    {
        std::cerr << "Error: Mismatched vector sizes." << std::endl;
        return false;
    }

    // Fill the map
    for (size_t i = 0; i < paramNames->size(); ++i)
    {
        branchMap[(*paramNames)[i]] = (*paramVals)[i];
        if (paramErr->size() > i)
            branchMap[(*paramNames)[i] + "_err"] = (*paramErr)[i];
    }

    branchMap["reactor_ratio"] = reacRatio->at(0);

    // Close the file
    file.Close();

    return true;
}

void makeFixedOscTree(const std::string &fitConfigFile_, const std::string &oscGridConfigFile_)
{
    Rand::SetSeed(0);

    // Load up the fit configuration information
    FitConfig fitConfig;
    FitConfigLoader fitLoader(fitConfigFile_);
    fitConfig = fitLoader.LoadActive();
    std::string outDir = fitConfig.GetOutDir();
    ParameterDict noms = fitConfig.GetNominals();
    ParameterDict constrMeans = fitConfig.GetConstrMeans();
    ParameterDict constrSigmas = fitConfig.GetConstrSigmas();
    ParameterDict mins = fitConfig.GetMinima();
    ParameterDict maxs = fitConfig.GetMaxima();

    for (ParameterDict::iterator it = noms.begin(); it != noms.end(); ++it)
    {
        if (!constrMeans[it->first])
        {
            constrMeans[it->first] = noms[it->first];
            constrSigmas[it->first] = -999;
        }
    }

    // Use the osc grid config to get the number of steps in each oscillation parameter
    OscGridConfigLoader oscGridLoader(oscGridConfigFile_);
    OscGridConfig oscGridConfig = oscGridLoader.Load();
    int numDeltam = oscGridConfig.GetNumValsDm21sq();
    int numTheta = oscGridConfig.GetNumValsSsqth12();

    double theta, deltam;

    // Create a ROOT file to store the tree
    std::string outFilename = "fit_result_tree.root";
    TFile outputFile((outDir + "/" + outFilename).c_str(), "RECREATE");
    TTree tree("fitResults", "Fit results from text files");

    // Make a map of branch names to values that will fill them
    std::map<std::string, double> branchMap = noms;

    std::vector<std::string> newKeys;
    for (const auto &pair : branchMap)
        newKeys.push_back(pair.first + "_err");
    // Now add the new entries
    for (const auto &key : newKeys)
        branchMap[key] = 0;

    branchMap["LLH"] = 0;
    branchMap["fit_valid"] = 0;
    branchMap["reactor_ratio"] = 0;

    // Initialise these to 0
    for (auto &[_, v] : branchMap)
        v = 0;

    // And tell it to fill with the values from the map
    for (auto &[name, value] : branchMap)
        tree.Branch(name.c_str(), &value);

    // Loop over theta and deltam
    for (int i = 0; i < numTheta; ++i)
    {
        theta = mins["theta12"] + i * ((maxs["theta12"] - mins["theta12"]) / numTheta);
        for (int j = 0; j < numDeltam; ++j)
        {
            deltam = mins["deltam21"] + j * ((maxs["deltam21"] - mins["deltam21"]) / numDeltam);

            // Construct file path
            std::ostringstream directory;
            directory << outDir << "/th" << std::fixed << std::setprecision(2) << theta
                      << "/th" << std::fixed << std::setprecision(2) << theta
                      << "_dm" << std::fixed << std::setprecision(8) << deltam;

            branchMap["theta12"] = theta;
            branchMap["deltam21"] = deltam;

            std::string filePath = directory.str() + "/fit_result.root";
            // Parse the file and fill the tree
            if (parseFitResultRoot(filePath, branchMap))
            {
                tree.Fill();
            }
        }
        std::cout << "done " << theta << std::endl;
    }

    // Write our tree and make vectors of the parameter names and asimov values
    outputFile.cd();
    tree.Write();
    std::vector<std::string> paramNameVec;
    std::vector<double> paramVals;
    branchMap = noms;
    mapToVectors(branchMap, paramNameVec, paramVals);
    outputFile.WriteObject(&paramNameVec, "param_names");
    outputFile.WriteObject(&paramVals, "param_asimov_values");

    std::vector<std::string> constrNameVec;
    std::vector<double> constrMeansVals;
    mapToVectors(constrMeans, constrNameVec, constrMeansVals);
    outputFile.WriteObject(&constrMeansVals, "constr_mean_values");

    std::vector<double> constrSigmaVals;
    mapToVectors(constrSigmas, constrNameVec, constrSigmaVals);
    outputFile.WriteObject(&constrSigmaVals, "constr_sigma_values");

    outputFile.Close();

    std::cout << "TTree saved to " << outDir << "/" << outFilename << std::endl;
}

int main(int argc, char *argv[])
{

    if (argc != 3)
    {
        std::cout << "\nUsage makeFixedOscTree <fit_config_file> <oscgrid_config_file>" << std::endl;
        return 1;
    }

    std::string fitConfigFile(argv[1]);
    std::string oscGridConfigFile(argv[2]);

    makeFixedOscTree(fitConfigFile, oscGridConfigFile);

    return 0;
}
