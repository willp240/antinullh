#include "Utilities.hh"

namespace antinufit
{

    // Taken from antinu rat-tools
    TVector3 LLAtoECEF(double longitude, double latitude, double altitude)
    {
        // reference http://www.mathworks.co.uk/help/aeroblks/llatoecefposition.html
        static double toRad = TMath::Pi() / 180.;
        static double Earthradius = 6378137.0; // Radius of the Earth (in meters)
        static double f = 1. / 298.257223563;  // Flattening factor WGS84 Model
        static double L, rs, x, y, z;
        L = atan(pow((1. - f), 2) * tan(latitude * toRad)) * 180. / TMath::Pi();
        rs = sqrt(pow(Earthradius, 2) / (1. + (1. / pow((1. - f), 2) - 1.) * pow(sin(L * toRad), 2)));
        x = (rs * cos(L * toRad) * cos(longitude * toRad) + altitude * cos(latitude * toRad) * cos(longitude * toRad)) / 1000; // in km
        y = (rs * cos(L * toRad) * sin(longitude * toRad) + altitude * cos(latitude * toRad) * sin(longitude * toRad)) / 1000; // in km
        z = (rs * sin(L * toRad) + altitude * sin(latitude * toRad)) / 1000;                                                   // in km

        TVector3 ECEF = TVector3(x, y, z);

        return ECEF;
    }

    // Taken from antinu rat-tools
    double GetReactorDistanceLLA(const double &longitude, const double &latitude, const double &altitude)
    {
        const TVector3 SNO_ECEF_coord_ = TVector3(672.87, -4347.18, 4600.51);
        double dist = (LLAtoECEF(longitude, latitude, altitude) - SNO_ECEF_coord_).Mag();
        return dist;
    }

    std::unordered_map<int, double> LoadIndexDistanceMap(std::string filename)
    {

        std::unordered_map<int, double> indexDistance;

        // Read the JSON file
        std::ifstream file(filename);
        if (!file.is_open())
        {
            std::cerr << "Utilities::LoadIndexDistanceMap() Could not open " << filename << std::endl;
            throw;
        }

        // Parse the JSON file into a json object
        nlohmann::json reactorData;
        file >> reactorData;

        // Loop through the JSON object and fill the maps
        for (const auto &[key, value] : reactorData.items())
        {
            int intKey = std::stoi(key); // Convert the reactor key to int
            double distance = value[1];  // Second element is the distance

            // Fill the maps
            indexDistance[intKey] = distance;
        }
        return indexDistance;
    }

    std::unordered_map<std::string, int> LoadNameIndexMap(std::string filename)
    {
        // First read the reactor distance info
        std::unordered_map<std::string, int> reactorNameIndex;

        // Read the JSON file
        std::ifstream file(filename);
        if (!file.is_open())
        {
            std::cerr << "Could not open reactors JSON file!" << std::endl;
            throw;
        }

        // Parse the JSON file into a json object
        nlohmann::json reactorData;
        file >> reactorData;

        // Loop through the JSON object and fill the maps
        for (const auto &[key, value] : reactorData.items())
        {
            int intKey = std::stoi(key);        // Convert the reactor key to int
            std::string reactorName = value[0]; // First element is the string

            // Fill the map
            reactorNameIndex[reactorName] = intKey;
        }

        return reactorNameIndex;
    }

    std::vector<double> LinSpace(double start, double end, size_t num_vals)
    {
        /*
         * Generates N linearly-spaced values between start and end.
         * From: https://stackoverflow.com/a/27030598
         */
        std::vector<double> linspaced;
        // = new std::vector<double>();

        const auto num = static_cast<double>(num_vals);

        if (num == 0)
        {
            return linspaced;
        }
        if (num == 1)
        {
            linspaced.push_back(start);
            return linspaced;
        }
        double delta = (end - start) / (num - 1.);

        for (size_t i = 0; i < num_vals - 1; ++i)
        {
            linspaced.push_back(start + delta * static_cast<double>(i));
        }
        linspaced.push_back(end); // I want to ensure that start and end
        // are exactly the same as the input
        return linspaced;
    }

    std::pair<size_t, size_t> GetLowerUpperIndices(const std::vector<double> vec, double val)
    {
        /*
         * Given an ordered vector of doubles, and a number lying between their min and max values,
         * return via a binary search the lower and upper indices of the vector that bound val.
         * E.g. if vec = {2, 3, 3.5, 4, 5}, val = 3.2, --> {1, 2}
         */
        const auto vec_high_it = std::lower_bound(vec.begin(), vec.end(), val);
        const size_t ind_high = std::distance(vec.begin(), vec_high_it);
        if (vec_high_it == vec.begin())
        {
            return {ind_high, ind_high + 1};
        }
        else
        {
            return {ind_high - 1, ind_high};
        }
    }

    std::vector<std::string> SplitString(const std::string &input, char delimiter)
    {
        std::vector<std::string> tokens;
        std::stringstream ss(input);
        std::string token;

        while (std::getline(ss, token, delimiter))
        {
            tokens.push_back(token);
        }
        return tokens;
    }

    void PrintParams(ParameterDict mins, ParameterDict maxs, ParameterDict noms, ParameterDict constrMeans, ParameterDict constrSigmas,
                     ParameterDict constrRatioMeans, ParameterDict constrRatioSigmas, std::map<std::string, std::string> constrRatioParName,
                     ParameterDict constrCorrs, std::map<std::string, std::string> constrCorrParName, std::map<std::string, std::vector<std::string>> datasets)
    {

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

        for (std::map<std::string, double>::iterator it = noms.begin(); it != noms.end(); it++)
        {
            if (std::find(tempNamesVec->begin(), tempNamesVec->end(), it->first) != tempNamesVec->end())
                continue;
            tempNamesVec->push_back(it->first);
        }

        std::cout << std::endl;
        std::cout << "************** Fit Parameters **************" << std::endl;
        std::cout << " -----------------------------------------------------------------------------------------------------------------------------------------------" << std::endl;
        std::cout << "| ";
        std::cout << std::left << std::setw(25) << "Name";
        std::cout << "| ";
        std::cout << std::left << std::setw(15) << "Nominal";
        std::cout << "| ";
        std::cout << std::left << std::setw(15) << "Minimum";
        std::cout << "| ";
        std::cout << std::left << std::setw(15) << "Maximum";
        std::cout << "| ";
        std::cout << std::left << std::setw(20) << "Datasets";
        std::cout << "| ";
        std::cout << std::left << std::setw(20) << "Constraint Mean";
        std::cout << "| ";
        std::cout << std::left << std::setw(20) << "Constraint Sigma";
        std::cout << "| " << std::endl;
        std::cout << " ===============================================================================================================================================" << std::endl;
        for (int iParam = 0; iParam < tempNamesVec->size(); iParam++)
        {

            if (!noms[tempNamesVec->at(iParam)])
                continue;

            std::cout << "| ";
            std::cout << std::left << std::setw(25) << tempNamesVec->at(iParam);
            std::cout << "| ";
            std::cout << std::left << std::setw(15) << noms[tempNamesVec->at(iParam)];
            std::cout << "| ";
            std::cout << std::left << std::setw(15) << mins[tempNamesVec->at(iParam)];
            std::cout << "| ";
            std::cout << std::left << std::setw(15) << maxs[tempNamesVec->at(iParam)];
            std::cout << "| ";
            std::string datasetString = "";
            for (int iDS = 0; iDS < datasets[tempNamesVec->at(iParam)].size(); iDS++)
            {
                datasetString += datasets[tempNamesVec->at(iParam)].at(iDS);
                if (iDS < datasets[tempNamesVec->at(iParam)].size() - 1)
                    datasetString += ", ";
            }
            std::cout << std::left << std::setw(20) << datasetString;
            std::cout << "| ";
            if (constrMeans[tempNamesVec->at(iParam)])
            {
                std::cout << std::left << std::setw(20) << constrMeans[tempNamesVec->at(iParam)];
                std::cout << "| ";
                std::cout << std::left << std::setw(20) << constrSigmas[tempNamesVec->at(iParam)];
                std::cout << "| ";
            }
            else
            {
                std::cout << std::left << std::setw(20) << "" << "| ";
                std::cout << std::left << std::setw(20) << "" << "| ";
            }
            std::cout << std::endl;
        }

        std::cout << " -----------------------------------------------------------------------------------------------------------------------------------------------" << std::endl;
        if (constrRatioMeans.size() > 0)
        {
            std::cout << std::endl;
            std::cout << "************** Ratio Constraints **************" << std::endl;
            std::cout << " -------------------------------------------------------------------------------------------------" << std::endl;
            std::cout << "| ";
            std::cout << std::left << std::setw(25) << "Parameter";
            std::cout << "| ";
            std::cout << std::left << std::setw(25) << "Parameter";
            std::cout << "| ";
            std::cout << std::left << std::setw(20) << "Constraint Mean";
            std::cout << "| ";
            std::cout << std::left << std::setw(20) << "Constraint Sigma";
            std::cout << "| " << std::endl;
            std::cout << " =================================================================================================" << std::endl;

            for (std::map<std::string, double>::iterator it = constrRatioMeans.begin(); it != constrRatioMeans.end(); it++)
            {
                std::cout << "| ";
                std::cout << std::left << std::setw(25) << it->first;
                std::cout << "| ";
                std::cout << std::left << std::setw(25) << constrRatioParName[it->first];
                std::cout << "| ";
                std::cout << std::left << std::setw(20) << constrRatioMeans[it->first];
                std::cout << "| ";
                std::cout << std::left << std::setw(20) << constrRatioSigmas[it->first];
                std::cout << "| " << std::endl;
            }
            std::cout << " -------------------------------------------------------------------------------------------------" << std::endl;
        }

        if (constrCorrs.size() > 0)
        {
            std::cout << std::endl;
            std::cout << "************** Correlations **************" << std::endl;
            std::cout << " ---------------------------------------------------------------------------" << std::endl;
            std::cout << "| ";
            std::cout << std::left << std::setw(25) << "Parameter";
            std::cout << "| ";
            std::cout << std::left << std::setw(25) << "Parameter";
            std::cout << "| ";
            std::cout << std::left << std::setw(20) << "Correlation";
            std::cout << "| ";
            std::cout << std::endl;
            std::cout << " ===========================================================================" << std::endl;

            for (std::map<std::string, double>::iterator it = constrCorrs.begin(); it != constrCorrs.end(); it++)
            {
                std::cout << "| ";
                std::cout << std::left << std::setw(25) << it->first;
                std::cout << "| ";
                std::cout << std::left << std::setw(25) << constrCorrParName[it->first];
                std::cout << "| ";
                std::cout << std::left << std::setw(20) << constrCorrs[it->first];
                std::cout << "| ";
            }
            std::cout << std::endl;
            std::cout << " ---------------------------------------------------------------------------" << std::endl;
        }

        std::cout << std::endl;
        return;
    }

    std::string stripQuoteMarks(std::string s_)
    {

        // Remove "" if it surrounds the tex name string
        size_t start = 0;
        size_t end = s_.size() - 1;
        if (s_[start] == '"' || s_[start] == '\'')
            start++;
        if (end > start && (s_[end] == '"' || s_[end] == '\''))
            end--;
        s_ = s_.substr(start, end - start + 1);

        return s_;
    }
}
