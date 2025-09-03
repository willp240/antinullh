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

// Function to overwrite the save_outputs bool and Minuit settings in a fit config file
void overwriteSaveOutputsAndMinuitSettings(const std::string &filepath)
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
    bool foundMinuitMethod = false;
    bool foundMinuitStrategy = false;
    bool foundMinuitTolerance = false;
    bool foundIterations = false;

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

        if (inSummary)
        {
            if (line.find("save_outputs") != std::string::npos)
            {
                line = "save_outputs = 1";
                foundSaveOutputs = true;
            }
            else if (line.find("minuit_method") != std::string::npos)
            {
                line = "minuit_method = \"Migrad\"";
                foundMinuitMethod = true;
            }
            else if (line.find("minuit_strategy") != std::string::npos)
            {
                line = "minuit_strategy = 2";
                foundMinuitStrategy = true;
            }
            else if (line.find("minuit_tolerance") != std::string::npos)
            {
                line = "minuit_tolerance = 0.01";
                foundMinuitTolerance = true;
            }
            else if (line.find("iterations") != std::string::npos && !foundIterations)
            {
                line = "iterations = 10000000";
                foundIterations = true;
            }
        }

        lines.push_back(line);
    }

    file.close();

    // Insert any missing entries into the [summary] section
    std::vector<std::string>::iterator lineIt = std::find(lines.begin(), lines.end(), "[summary]");
    if (lineIt != lines.end())
    {
        ++lineIt; // move to line after [summary]

        if (!foundSaveOutputs)
            lineIt = lines.insert(lineIt, "save_outputs = 1") + 1;
        if (!foundMinuitMethod)
            lineIt = lines.insert(lineIt, "minuit_method = \"Migrad\"") + 1;
        if (!foundMinuitStrategy)
            lineIt = lines.insert(lineIt, "minuit_strategy = 2") + 1;
        if (!foundMinuitTolerance)
            lineIt = lines.insert(lineIt, "minuit_tolerance = 0.01") + 1;
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
    std::cout << "Updated summary section in " << filepath << std::endl;
}


// Function to parse output file from job of many fits to a map
bool parseFitResultsTxt(const std::string &filename, std::map<std::string, double> &branchMap,
                        TTree *tree, double *bestLLH, double *bestDeltam, double *bestTheta, std::string theta12name)
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
        if (line.find(theta12name.c_str()) != std::string::npos)
        {
            std::istringstream iss(line);
            std::string dummy;
            double theta;
            iss >> dummy >> theta;
            branchMap[theta12name] = theta;
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
                *bestTheta = branchMap[theta12name];
            }
        }
        if (line.find("FitValid:") != std::string::npos)
        {
            std::istringstream iss(line);
            std::string dummy;
            double fitValid;
            iss >> dummy >> fitValid;
            branchMap["fit_valid"] = fitValid;

            // FitValid is last line for this fit, so now let's fill
            tree->Fill();
        }
        if (line.find("Reactor Ratio:") != std::string::npos)
        {

            std::istringstream iss(line);
            std::string datasetName, dummy;
            double reactorRatio;

            // Format: datasetname Reactor Ratio: value
            iss >> datasetName >> dummy >> dummy; // Skip "Reactor" and "Ratio:"
            iss >> reactorRatio;

            std::string key = datasetName + "_ratio";
            branchMap[key] = reactorRatio;
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

    for (ParameterDict::iterator parIt = noms.begin(); parIt != noms.end(); ++parIt)
    {
        if (!constrMeans[parIt->first])
        {
            constrMeans[parIt->first] = noms[parIt->first];
            constrSigmas[parIt->first] = -999;
        }
    }

    // Check we've got the oscillation parameters we need, and whether theta is sin-ed or not
    bool hasTheta12 = noms.find("theta12") != noms.end();
    bool hasSinTheta12 = noms.find("sintheta12") != noms.end();
    bool hasSinSqTheta12 = noms.find("sinsqtheta12") != noms.end();
    bool hasDeltam21 = noms.find("deltam21") != noms.end();
    int thetacount = int(hasTheta12) + int(hasSinTheta12) + int(hasSinSqTheta12);
    if (thetacount == 0 || !hasDeltam21)
    {
        throw std::runtime_error("ERROR: A theta12, sintheta12, or sinsqtheta12 parameter, along with a deltam21 parameter, must be provided in the fitconfig.");
    }
    if (thetacount > 1)
    {
        throw std::runtime_error("ERROR: More than one of theta12, sintheta12, sinsqtheta12 parameters were set in the fitconfig. There must be exactly one.");
    }
    std::string theta12name;
    if (hasTheta12)
        theta12name = "theta12";
    else if (hasSinTheta12)
        theta12name = "sintheta12";
    else if (hasSinSqTheta12)
        theta12name = "sinsqtheta12";

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

    // Make error entry for each parameter
    std::vector<std::string> newKeys;
    for (const auto &pair : branchMap)
    {
        newKeys.push_back(pair.first + "_err");
        if (pair.first.find("reactor") != std::string::npos)
        {
            newKeys.push_back(pair.first + "_ratio");
        }
    }
    // Now add the new entries
    for (const auto &key : newKeys)
        branchMap[key] = 0;

    branchMap["LLH"] = 0;
    branchMap["fit_valid"] = 0;

    // Make sure everything's Initialised to 0
    for (auto &[_, v] : branchMap)
        v = 0;

    // And tell it to fill with the values from the map
    for (auto &[name, value] : branchMap)
        tree.Branch(name.c_str(), &value);

    double bestLLH = 999;
    double bestDeltam = 999;
    double bestTheta = 999;

    // Loop over all .output files in output dir
    for (const auto &entry : std::filesystem::directory_iterator((outDir + "/output").c_str()))
    {
        if (entry.is_regular_file() && entry.path().extension() == ".output")
        {
            std::string filePath = entry.path().string();
            std::cout << "Parsing " << filePath << std::endl;

            parseFitResultsTxt(filePath, branchMap, &tree, &bestLLH, &bestDeltam, &bestTheta, theta12name);
        }
    }

    // Write our tree and make vectors of the parameter names and asimov values
    outputFile.cd();
    tree.Write();


    std::cout << std::endl << "TTree saved to " << outDir << "/" << outFilename << std::endl << std::endl;

    std::cout << "Now rerunning fit with save outputs flag on for deltam: " << bestDeltam << ", " << theta12name << ": " << bestTheta << std::endl
              << std::endl;

    // Find fit cfg path for that fit
    std::ostringstream bestfitcfg;
    bestfitcfg << outDir << "/th" << std::fixed << std::setprecision(3) << bestTheta << "/cfg/" << "fit_config_"
               << "th" << std::fixed << std::setprecision(3) << bestTheta
               << "_dm" << std::fixed << std::setprecision(8) << bestDeltam << ".ini";

    // Overwrite save_output bool and Minuit settings in that config
    overwriteSaveOutputsAndMinuitSettings(bestfitcfg.str());

    std::ostringstream eventcfg;
    eventcfg << outDir << "/th" << std::fixed << std::setprecision(3) << bestTheta << "/cfg/" << "event_config.ini";
    std::ostringstream pdfcfg;
    pdfcfg << outDir << "/th" << std::fixed << std::setprecision(3) << bestTheta << "/cfg/" << "pdf_config.ini";
    std::ostringstream syscfg;
    syscfg << outDir << "/th" << std::fixed << std::setprecision(3) << bestTheta << "/cfg/" << "syst_config.ini";
    std::ostringstream osccfg;
    osccfg << outDir << "/th" << std::fixed << std::setprecision(3) << bestTheta << "/cfg/" << "oscgrid_config.ini";

    std::ostringstream fit_command;
    fit_command << "./bin/fixedosc_fit " << bestfitcfg.str() << " " << eventcfg.str() << " " << pdfcfg.str() << " " << syscfg.str() << " " << osccfg.str();
    std::cout << fit_command.str().c_str() << std::endl;
    // Now rerun that fit
    system(fit_command.str().c_str());

    // Now we open that fit file, get the par name vec, par vals vec, and par err vec, and cov matrix
    // These will all have the parameters in the same order
    std::ostringstream fitfilename;
    fitfilename << outDir << "/th" << std::fixed << std::setprecision(3) << bestTheta << "/th" << std::fixed << std::setprecision(3)
                << bestTheta << "_dm" << std::fixed << std::setprecision(8) << bestDeltam << "/fit_result.root";

    // Open the ROOT file
    TFile *fitfile = TFile::Open(fitfilename.str().c_str(), "READ");
    if (!fitfile || fitfile->IsZombie())
    {
        std::cerr << "Error: Could not open file " << fitfilename.str() << std::endl;
        return;
    }

    // Get the covariance TMatrixT
    TMatrixT<double> *covMatrix = nullptr;
    fitfile->GetObject("covMatrix", covMatrix);
    if (!covMatrix)
    {
        std::cerr << "Error: Could not find TMatrixT 'covMatrix' in file " << fitfilename.str() << std::endl;
        fitfile->Close();
        return;
    }

    // Get the parameter name vector
    std::vector<std::string> *paramNameVec = nullptr;
    fitfile->GetObject("paramNames", paramNameVec);
    if (!paramNameVec)
    {
        std::cerr << "Error: Could not find 'paramNames' in file " << fitfilename.str() << std::endl;
        fitfile->Close();
        return;
    }
    // Remove LLH and FitValid elements
    paramNameVec->resize(paramNameVec->size()-2);

    // Get the parameter values vector
    std::vector<double> *paramVals = nullptr;
    fitfile->GetObject("paramVals", paramVals);
    if (!paramVals)
    {
        std::cerr << "Error: Could not find 'paramVals' in file " << fitfilename.str() << std::endl;
        fitfile->Close();
        return;
    }
    // Remove LLH and FitValid elements
    paramVals->resize(paramVals->size()-2);


    // Get the parameter errors vector
    std::vector<double> *paramErrs = nullptr;
    fitfile->GetObject("paramErr", paramErrs);
    if (!paramErrs)
    {
        std::cerr << "Error: Could not find 'paramErrs' in file " << fitfilename.str() << std::endl;
        fitfile->Close();
        return;
    }

    // Now we're going to loop over the paramNameVec, and get the nominals, constraints and labels and put them in a vector
    // so it's in the same order as the cov matrix, paramVals, and paramErrs
    std::vector<double> nomVals;
    std::vector<double> constrMeansVals;
    std::vector<double> constrSigmaVals;
    std::vector<std::string> labelsVec;
    branchMap = noms;

    // Loop over paramNames
    for (int iPar = 0; iPar < paramNameVec->size(); iPar++)
    {
        nomVals.push_back(noms[paramNameVec->at(iPar)]);
        labelsVec.push_back(labels[paramNameVec->at(iPar)]);
        constrMeansVals.push_back(constrMeans[paramNameVec->at(iPar)]);
        constrSigmaVals.push_back(constrSigmas[paramNameVec->at(iPar)]);
    }

    // Now we have all the vectors are in the same order as the covariance matrix, 
    // we'll add the osc pars to the end of the vector
    paramNameVec->push_back("deltam21");
    paramNameVec->push_back(theta12name);
    paramVals->push_back(bestDeltam);
    paramVals->push_back(bestTheta);
    nomVals.push_back(noms["deltam21"]);
    nomVals.push_back(noms[theta12name]);
    labelsVec.push_back(labels["deltam21"]);
    labelsVec.push_back(labels[theta12name]);
    constrMeansVals.push_back(constrMeans["deltam21"]);
    constrMeansVals.push_back(constrMeans[theta12name]);
    constrSigmaVals.push_back(constrSigmas["deltam21"]);
    constrSigmaVals.push_back(constrSigmas[theta12name]);

    outputFile.WriteObject(paramNameVec, "param_names");
    outputFile.WriteObject(paramVals, "param_fit_values");
    outputFile.WriteObject(paramErrs, "param_fit_err");
    outputFile.WriteObject(&nomVals, "param_nom_values");
    outputFile.WriteObject(&constrMeansVals, "constr_mean_values");
    outputFile.WriteObject(&constrSigmaVals, "constr_sigma_values");
    outputFile.WriteObject(&labelsVec, "tex_labels");

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
