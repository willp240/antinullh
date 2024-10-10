#include "Utilities.hh"

namespace antinufit
{

    std::vector<double> linspace(double start, double end, size_t num_vals)
    {
        /*
         * Generates N linearly-spaced values between start and end.
         * From: https://stackoverflow.com/a/27030598
         */
        std::vector<double> linspaced;
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
    }

}
