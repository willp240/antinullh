#include "PlotFitParams.hh"
#include <stdexcept>

void plotMinuitParams(const char *filename)
{
    TFile file(filename, "READ");
    if (file.IsZombie())
        throw std::runtime_error("Cannot open Minuit postfit metadata");
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
        throw std::runtime_error("Incomplete Minuit postfit metadata");
    const size_t n = names->size();
    if (!n || values->size() != n || errors->size() != n || nominals->size() != n ||
        means->size() != n || sigmas->size() != n || labels->size() != n ||
        covariance->GetNrows() != n || covariance->GetNcols() != n)
        throw std::runtime_error("Inconsistent Minuit postfit metadata dimensions");
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
                               means, sigmas, covariance, theta, correlated, classifier);
}
