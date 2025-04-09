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
#include <filesystem>

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

// Function to overwrite the save_outputs bool in a fit config file
void overwriteSaveOutputsBool(const std::string &filepath)
{
    std::ifstream file(filepath);
    if (!file)
    {
        std::cerr << "Error: Could not open file " << filepath << std::endl;
        return;
    }
    std::vector<std::string> lines;
    std::string line;
    bool inSummary = false;
    bool foundSaveOutputs = false;

    while (std::getline(file, line))
    {
        // Check if we are in the [summary] section
        if (line == "[summary]")
        {
            inSummary = true;
        }
        else if (line.find("[") == 0 && line != "[summary]")
        {
            inSummary = false;
        }

        // Modify save_outputs if found in [summary]
        if (inSummary && line.find("save_outputs") != std::string::npos)
        {
            line = "save_outputs = 1";
            foundSaveOutputs = true;
        }

        lines.push_back(line);
    }

    file.close();

    // If save_outputs was not found, add it to the summary section
    if (!foundSaveOutputs)
    {
        std::vector<std::string>::iterator it = std::find(lines.begin(), lines.end(), "[summary]");
        it++;
        lines.insert(it, "save_outputs = 1");
    }

    // Write back to file
    std::ofstream outFile(filepath);
    if (!outFile)
    {
        std::cerr << "Error: Could not write to file " << filepath << std::endl;
        return;
    }

    for (const std::string &l : lines)
    {
        outFile << l << "\n";
    }

    outFile.close();
    std::cout << "Updated save_outputs in " << filepath << std::endl;
}

// Function to turn a map of string to double into two vectors, one of strings, one of doubles
// This is useful for writing to root files
void doubleMapToVectors(std::map<std::string, double> sdmap, std::vector<std::string> &keysvec, std::vector<double> &valsvec)
{

    for (std::map<std::string, double>::iterator it = sdmap.begin(); it != sdmap.end(); ++it)
    {
        keysvec.push_back(it->first);
        valsvec.push_back(it->second);
    }
}

// Function to turn a map of string to string into two vectors of strings
// This is useful for writing to root files
void stringMapToVectors(std::map<std::string, std::string> sdmap, std::vector<std::string> &keysvec, std::vector<std::string> &valsvec)
{

    for (std::map<std::string, std::string>::iterator it = sdmap.begin(); it != sdmap.end(); ++it)
    {
        keysvec.push_back(it->first);
        valsvec.push_back(it->second);
    }
}

// Function to parse output file from job of many fits to a map
bool parseFitResultsTxt(const std::string &filename, std::map<std::string, double> &branchMap,
                        TTree *tree, double *bestLLH, double *bestDeltam, double *bestTheta)
{

    std::ifstream infile(filename);
    if (!infile)
    {
        std::cerr << "Failed to open file: " << filename << std::endl;
        return false;
    }

    bool inParamsSection = false;

    std::string line;
    while (std::getline(infile, line))
    {
        if (line.find("Best Fit Values:") != std::string::npos)
        {
            inParamsSection = true;
            continue;
        }

        // If we reach an empty line, we exit the parameters section
        if (inParamsSection && line.find("Fit complete for:") != std::string::npos)
        {
            inParamsSection = false;
            continue;
        }

        // Extract parameter names and values dynamically
        if (inParamsSection)
        {
            std::istringstream iss(line);
            std::string paramName;
            double value;

            // Read entire line until the last token, which should be the parameter value
            std::string word;
            while (iss >> word)
            {
                try
                {
                    // Try converting the last token to a double
                    value = std::stod(word);
                    break;
                }
                catch (const std::invalid_argument &)
                {
                    // Not a number, append to parameter name
                    if (!paramName.empty())
                        paramName += " ";
                    paramName += word;
                }
            }
            if (!paramName.empty())
            {
                branchMap[paramName] = value;
            }
        }

        // Explicitly read "deltam", "theta", "LLH", and "Fit Valid"
        if (line.find("deltam") != std::string::npos)
        {
            std::istringstream iss(line);
            std::string dummy;
            double deltam;
            iss >> dummy >> deltam;
            branchMap["deltam21"] = deltam;
        }
        if (line.find("theta") != std::string::npos)
        {
            std::istringstream iss(line);
            std::string dummy;
            double theta;
            iss >> dummy >> theta;
            branchMap["theta12"] = theta;
        }
        if (line.find("LLH:") != std::string::npos)
        {
            std::istringstream iss(line);
            std::string dummy;
            double llh;
            iss >> dummy >> llh;
            branchMap["LLH"] = llh;
            if (branchMap["LLH"] < *bestLLH)
            {
                *bestLLH = branchMap["LLH"];
                *bestDeltam = branchMap["deltam21"];
                *bestTheta = branchMap["theta12"];
            }
        }
        if (line.find("FitValid:") != std::string::npos)
        {
            std::istringstream iss(line);
            std::string dummy;
            double fitValid;
            iss >> dummy >> fitValid;
            branchMap["fit_valid"] = fitValid;
        }
        if (line.find("ReactorRatio:") != std::string::npos)
        {
            std::istringstream iss(line);
            std::string dummy;
            double reactorRatio;
            iss >> dummy >> reactorRatio;
            branchMap["reactor_ratio"] = reactorRatio;
            // ReactorRatio is last line for this fit, so now let's fill
            tree->Fill();
        }
    }
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
    std::map<std::string, std::string> labels = fitConfig.GetTexLabels();

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

    double bestLLH = 999;
    double bestDeltam = 999;
    double bestTheta = 999;

    // Loop over theta and deltam
    for (int i = 0; i < numTheta; ++i)
    {
        theta = mins["theta12"] + i * ((maxs["theta12"] - mins["theta12"]) / numTheta);

        // Construct file path
        std::string dirName = std::filesystem::path(outDir).filename().string();
        std::ostringstream filePath;
        filePath << outDir << "/output/" << dirName
                 << "_th" << std::fixed << std::setprecision(2) << theta
                 << ".output";

        // Parse the file and fill the tree
        parseFitResultsTxt(filePath.str(), branchMap, &tree, &bestLLH, &bestDeltam, &bestTheta);
        std::cout << "Done " << theta << std::endl;
    }

    // Write our tree and make vectors of the parameter names and asimov values
    outputFile.cd();
    tree.Write();
    std::vector<std::string> paramNameVec;
    std::vector<double> paramVals;
    branchMap = noms;
    doubleMapToVectors(branchMap, paramNameVec, paramVals);
    outputFile.WriteObject(&paramNameVec, "param_names");
    outputFile.WriteObject(&paramVals, "param_asimov_values");

    std::vector<std::string> constrNameVec;
    std::vector<double> constrMeansVals;
    doubleMapToVectors(constrMeans, constrNameVec, constrMeansVals);
    outputFile.WriteObject(&constrMeansVals, "constr_mean_values");

    std::vector<double> constrSigmaVals;
    doubleMapToVectors(constrSigmas, constrNameVec, constrSigmaVals);
    outputFile.WriteObject(&constrSigmaVals, "constr_sigma_values");

    std::vector<std::string> labelsVals;
    stringMapToVectors(labels, paramNameVec, labelsVals);
    outputFile.WriteObject(&labelsVals, "tex_labels");

    std::cout << "TTree saved to " << outDir << "/" << outFilename << std::endl;

    std::cout << "Now rerunning fit with save outputs flag on for deltam: " << bestDeltam << ", theta: " << bestTheta << std::endl
              << std::endl;

    // Find fit cfg path for that fit
    std::ostringstream bestfitcfg;
    bestfitcfg << outDir << "/th" << std::fixed << std::setprecision(2) << bestTheta << "/cfg/" << "fit_config_"
               << "th" << std::fixed << std::setprecision(2) << bestTheta
               << "_dm" << std::fixed << std::setprecision(8) << bestDeltam << ".ini";

    // Overwrite save_output bool in that config
    overwriteSaveOutputsBool(bestfitcfg.str());

    std::ostringstream eventcfg;
    eventcfg << outDir << "/th" << std::fixed << std::setprecision(2) << bestTheta << "/cfg/" << "event_config.ini";
    std::ostringstream pdfcfg;
    pdfcfg << outDir << "/th" << std::fixed << std::setprecision(2) << bestTheta << "/cfg/" << "pdf_config.ini";
    std::ostringstream syscfg;
    syscfg << outDir << "/th" << std::fixed << std::setprecision(2) << bestTheta << "/cfg/" << "syst_config.ini";
    std::ostringstream osccfg;
    osccfg << outDir << "/th" << std::fixed << std::setprecision(2) << bestTheta << "/cfg/" << "oscgrid_config.ini";

    std::ostringstream fit_command;
    fit_command << "./bin/fixedosc_fit " << bestfitcfg.str() << " " << eventcfg.str() << " " << pdfcfg.str() << " " << syscfg.str() << " " << osccfg.str();

    // Now rerun that fit
    system(fit_command.str().c_str());

    std::ostringstream fitfilename;
    fitfilename << outDir << "/th" << std::fixed << std::setprecision(2) << bestTheta << "/th" << std::fixed << std::setprecision(2)
                << bestTheta << "_dm" << std::fixed << std::setprecision(8) << bestDeltam << "/fit_result.root";

    // Open the ROOT file
    TFile *fitfile = TFile::Open(fitfilename.str().c_str(), "READ");
    if (!fitfile || fitfile->IsZombie())
    {
        std::cerr << "Error: Could not open file " << fitfilename.str() << std::endl;
        return;
    }

    // Get the TMatrixT
    TMatrixT<double> *covMatrix = nullptr;
    fitfile->GetObject("covMatrix", covMatrix);
    if (!covMatrix)
    {
        std::cerr << "Error: Could not find TMatrixT 'covMatrix' in file " << fitfilename.str() << std::endl;
        fitfile->Close();
        return;
    }

    outputFile.cd();
    covMatrix->Write("covMatrix");
    outputFile.Close();
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
