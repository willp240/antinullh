#include "FitParamDrawing.hh"
#include <stdexcept>

void plotFitParams(const char *filename, bool posterior = false)
{
    TFile file(filename, "READ");
    if (file.IsZombie())
        throw std::runtime_error("Cannot open postfit summary metadata");
    std::vector<std::string> *names = nullptr;
    std::vector<std::string> *labels = nullptr;
    std::vector<double> *values = nullptr;
    std::vector<double> *errors = nullptr;
    std::vector<double> *nominals = nullptr;
    std::vector<double> *means = nullptr;
    std::vector<double> *sigmas = nullptr;
    TMatrixT<double> *covariance = nullptr;
    file.GetObject("param_names", names);
    file.GetObject("tex_labels", labels);
    file.GetObject("param_fit_values", values);
    file.GetObject("param_fit_err", errors);
    file.GetObject("param_nom_values", nominals);
    file.GetObject("constr_mean_values", means);
    file.GetObject("constr_sigma_values", sigmas);
    file.GetObject("covMatrix", covariance);
    if (!names || !labels || !values || !errors || !nominals || !means || !sigmas || !covariance)
        throw std::runtime_error("Incomplete postfit summary metadata");
    const size_t n = names->size();
    if (!n || values->size() != n || errors->size() != n || nominals->size() != n ||
        means->size() != n || sigmas->size() != n || labels->size() != n ||
        covariance->GetNrows() != n || covariance->GetNcols() != n)
        throw std::runtime_error("Inconsistent postfit summary metadata dimensions");
    std::map<std::string, std::pair<double, double>> intervals;
    std::vector<antinupostfit::ParameterEstimate> estimates;
    bool modePosterior = false;
    if (posterior)
    {
        std::vector<double> *modes = nullptr;
        file.GetObject("posterior_mode", modes);
        if (modes)
        {
            if (modes->size() != n)
                throw std::runtime_error("Inconsistent posterior mode dimensions");
            std::vector<std::string> *regionNames = nullptr;
            std::vector<double> *regionLow = nullptr;
            std::vector<double> *regionHigh = nullptr;
            file.GetObject("posterior_mode_hpd_names", regionNames);
            file.GetObject("posterior_mode_hpd_low", regionLow);
            file.GetObject("posterior_mode_hpd_high", regionHigh);
            if (!regionNames || !regionLow || !regionHigh ||
                regionNames->size() != regionLow->size() || regionNames->size() != regionHigh->size())
                throw std::runtime_error("Missing or inconsistent highest-density region metadata");
            antinupostfit::ParameterEstimate estimate = {"posterior_mode", "Mode, 68% Highest Density", kRed, 20, {}, {}};
            for (size_t i = 0; i < regionNames->size(); ++i)
            {
                if (!std::isfinite((*regionLow)[i]) || !std::isfinite((*regionHigh)[i]) || (*regionLow)[i] > (*regionHigh)[i])
                    throw std::runtime_error("Invalid highest-density region");
                estimate.regions[(*regionNames)[i]].push_back(std::make_pair((*regionLow)[i], (*regionHigh)[i]));
            }
            for (size_t i = 0; i < n; ++i)
            {
                bool containsMode = false;
                for (const std::pair<double, double> &region : estimate.regions[(*names)[i]])
                    if (region.first <= (*modes)[i] && (*modes)[i] <= region.second)
                    {
                        containsMode = true;
                        intervals[(*names)[i]] = std::make_pair((*modes)[i] - region.first, region.second - (*modes)[i]);
                    }
                if (!std::isfinite((*modes)[i]) || !containsMode)
                    throw std::runtime_error("Mode is missing from its highest-density regions");
                estimate.values[(*names)[i]] = std::make_pair((*modes)[i], 0.0);
            }
            values = modes;
            modePosterior = true;
            estimates.push_back(estimate);
            const std::vector<std::string> keys = {"param_fit_values", "posterior_mean", "posterior_gaussian_mean", "posterior_best"};
            const std::vector<std::string> widthKeys = {"", "posterior_sd", "posterior_gaussian_sigma", ""};
            const std::vector<std::string> plotKeys = {"posterior_median", "posterior_mean", "posterior_gaussian", "posterior_best"};
            const std::vector<std::string> descriptions = {"Median, 68% Equal-Tailed", "Mean, Standard Deviation", "Gaussian Mean, Fitted Width", "Maximum LLH"};
            const std::vector<int> colors = {kRed, kOrange + 7, kGreen + 2, kMagenta + 1};
            const std::vector<int> markers = {24, 21, 22, 33};
            std::vector<double> *medianLow = nullptr;
            std::vector<double> *medianHigh = nullptr;
            file.GetObject("posterior_error_low", medianLow);
            file.GetObject("posterior_error_high", medianHigh);
            if (!medianLow || !medianHigh || medianLow->size() != n || medianHigh->size() != n)
                throw std::runtime_error("Missing median interval metadata");
            for (size_t j = 0; j < keys.size(); ++j)
            {
                std::vector<double> *central = nullptr;
                std::vector<double> *width = nullptr;
                file.GetObject(keys[j].c_str(), central);
                if (!widthKeys[j].empty())
                    file.GetObject(widthKeys[j].c_str(), width);
                if (!central)
                    continue;
                if (central->size() != n || (!widthKeys[j].empty() && (!width || width->size() != n)))
                    throw std::runtime_error("Inconsistent posterior estimate dimensions");
                antinupostfit::ParameterEstimate alternative = {plotKeys[j], descriptions[j], colors[j], markers[j], {}, {}};
                for (size_t i = 0; i < n; ++i)
                {
                    alternative.values[(*names)[i]] = std::make_pair((*central)[i], width ? (*width)[i] : 0.0);
                    if (plotKeys[j] == "posterior_median")
                        alternative.regions[(*names)[i]].push_back(std::make_pair((*central)[i] - (*medianLow)[i], (*central)[i] + (*medianHigh)[i]));
                }
                estimates.push_back(alternative);
            }
        }
        else
        {
            // Older files without modes retain their saved median intervals
            std::vector<double> *low = nullptr;
            std::vector<double> *high = nullptr;
            file.GetObject("posterior_error_low", low);
            file.GetObject("posterior_error_high", high);
            if (!low || !high || low->size() != n || high->size() != n)
                throw std::runtime_error("Missing posterior interval metadata");
            for (size_t i = 0; i < n; ++i)
                intervals[(*names)[i]] = std::make_pair((*low)[i], (*high)[i]);
        }
    }
    std::map<std::string, double> parameters;
    bool correlated = false;
    bool classifier = false;
    std::string theta;
    for (size_t i = 0; i < n; ++i)
    {
        parameters[(*names)[i]] = (*values)[i];
        if ((*names)[i] == "reactor_nubar_norm")
            correlated = true;
        if ((*names)[i].find("class_") == 0)
            classifier = true;
        if ((*names)[i] == "theta12" || (*names)[i] == "sintheta12" || (*names)[i] == "sinsqtheta12")
            theta = (*names)[i];
    }
    // Rates in a floated model are unoscillated before applying the shape systematic
    parameters["reactor_nubar_ratio"] = 1;
    parameters["reactor_nubar2_ratio"] = 1;
    antinupostfit::DrawFitParams(filename, parameters, names, errors, labels, nominals,
                               means, sigmas, covariance, theta, correlated, classifier, posterior ? &intervals : nullptr,
                               estimates.empty() ? nullptr : &estimates, modePosterior);
}
