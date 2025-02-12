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

// Function to turn a map of string to double into two vectors, one of strings one of doubles.
// The element number of one vector should correspond to the element number in the other.
// This is useful for writing to root files
void mapToVectors(std::map<std::string, double> sdmap, std::vector<std::string> &stringvec, std::vector<double> &doublevec)
{

    for (std::map<std::string, double>::iterator it = sdmap.begin(); it != sdmap.end(); ++it)
    {
        stringvec.push_back(it->first);
        doublevec.push_back(it->second);
    }
}

// Function to parse an outputted txt file and extract the fit parameters
bool ParseFitResult(const std::string &filename, std::map<std::string, double> &branchMap)
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

        // TO DO loop over params here
        if (iss >> key >> value)
        {
            branchMap[key] = value;
            if (key == "LLH:")
                branchMap["LLH"] = value;
        }
    }

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
    ParameterDict mins = fitConfig.GetMinima();
    ParameterDict maxs = fitConfig.GetMaxima();

    // Use the osc grid config to get the number of steps in each oscillation parameter
    OscGridConfigLoader oscGridLoader(oscGridConfigFile_);
    OscGridConfig oscGridConfig = oscGridLoader.Load();
    int numDeltam = oscGridConfig.GetNumValsDm21sq();
    int numTheta = oscGridConfig.GetNumValsSsqth12();

    double theta, deltam;

    // Create a ROOT file to store the tree
    std::string outFilename = "fit_results.root";
    TFile outputFile((outDir + "/" + outFilename).c_str(), "RECREATE");
    TTree tree("FitResults", "Fit results from text files");

    // Make a map of branch names to values that will fill them
    std::map<std::string, double> branchMap = noms;
    branchMap["LLH"] = 0;

    // Initialise these to 0
    for (auto &[_, v] : branchMap)
        v = 0;
    
    // And tell it to fill with the values from the map
    for (auto &[name, value] : branchMap)
        tree.Branch(name.c_str(), &value);

    // Oscillation parameters aren't read from the output file
    tree.Branch("theta", &theta);
    tree.Branch("deltam", &deltam);

    // Loop over theta and deltam
    // Read max and mins from config
    for (int i = 0; i < numTheta; ++i)
    {
        theta = mins["theta12"] + i * ((maxs["theta12"] - mins["theta12"]) / numTheta);
        for (int j = 0; j < numDeltam; ++j)
        {
            deltam = mins["deltam21"] + j * ((maxs["deltam21"] - mins["deltam21"]) / numDeltam);

            // Construct file path
            std::ostringstream directory;
            directory << outDir << "/th" << std::fixed << std::setprecision(2) << theta
                      << "_dm" << std::fixed << std::setprecision(8) << deltam;
            std::string filePath = directory.str() + "/fit_result.txt";
           
            // Parse the file and fill the tree
            if (ParseFitResult(filePath, branchMap))
            {
                tree.Fill();
            }
        }
        std::cout << "done " << theta << std::endl;
    }

    // Write our tree and make vectors of the parameter names and asimov values
    tree.Write();
    std::vector<std::string> paramNameVec;
    std::vector<double> paramVals;
    mapToVectors(branchMap, paramNameVec, paramVals);
    outputFile.WriteObject(&paramNameVec, "param_names");
    outputFile.WriteObject(&paramVals, "param_asimov_values");

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
