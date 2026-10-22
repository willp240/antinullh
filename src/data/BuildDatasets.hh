#pragma once

#include <map>
#include <string>
#include <vector>
#include <unordered_map>

#include <DatasetSetup.hh>
#include <BuildSysts.hh>
#include <OscParams.hh>
#include <FitParameters.hh>

#include <BinnedNLLH.h>
#include <BinnedED.h>

namespace antinufit
{
    // Match a component distribution to the observables used by the dataset
    inline BinnedED ProjectToDataObs(const BinnedED &src,
                                     const BinnedED &dataDist,
                                     const std::vector<std::string> &dataObs)
    {
        if (src.GetNDims() != dataDist.GetNDims())
            return src.Marginalise(dataObs);

        return src;
    }

    // Likelihood inputs, selected data, and optional distributions for one dataset
    struct PreparedDataset
    {
        // Inputs to the per-dataset likelihood
        std::vector<BinnedED> pdfs;
        std::vector<std::vector<std::string>> pdfGroups;
        std::vector<int> genRates;
        std::vector<NormFittingStatus> normFittingStatuses;
        BinnedED dataDist;
        ParameterDict initialValues;

        // Nominal and fake-data counts used to select or write the data distribution
        BinnedED asimov;
        BinnedED fakeDataset;
        std::vector<BinnedED> nominalComponents; // Retained only for requested output
        std::vector<BinnedED> fakeComponents;    // Retained only for requested output
        std::map<std::string, double> nominalNonBufferRates;

        // Reactor factors calculated while filling the component histograms
        std::map<std::string, double> reactorRatio;
        std::map<std::string, double> reactorFakeRatio;
    };

    // Lookup from configured dataset name to its prepared data
    using PreparedDatasetMap = std::map<std::string, PreparedDataset>;

    // Scale fit rates using the reactor ratios measured for all datasets
    void ApplyReactorRateScales(FitParameters &p,
                                const PreparedDatasetMap &datasets);

    PreparedDatasetMap BuildDatasets(
        const FitConfig &fitConfig, const FitParameters &p,
        const DatasetSetup &setup,
        const SystsByDataset &syst,
        const SystConfig &systConfig,
        const OscParams &osc,
        const std::unordered_map<int, double> &indexDistance,
        const DatasetParameterMap &datasetPars);

} // namespace antinufit
