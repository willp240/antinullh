#pragma once

#include <map>
#include <string>
#include <vector>
#include <unordered_map>

#include <FitInputs.hh>
#include <DatasetSetup.hh>
#include <SystSetup.hh>
#include <OscParams.hh>
#include <OutputDirs.hh>
#include <DatasetModel.hh>
#include <ParamUtils.hh>
#include <BuildComponent.hh>

#include <BinnedNLLH.h>
#include <BinnedED.h>

namespace antinufit
{
    struct BuildDatasetsConfig
    {
        BuildComponentConfig component;
    };

    // All the info needed to bring PDFs, systs, and datasets together to make the llh
    struct BuildResult
    {
        DatasetModelMap models;
        std::vector<BinnedNLLH> llhs;
        std::map<std::string, ParameterDict> initialByDataset;
        std::map<std::string, BinnedED> dataDists;
        std::map<std::string, double> reactorRatio;
        std::map<std::string, double> reactorRatioFD;
        std::map<std::string, std::map<std::string, double>> nominalNonBufferRates;

        // Optional scan payload for future llh-scan style executables.
        // dataset -> pdf -> [one entry per configured scan point]
        std::map<std::string, std::map<std::string, std::vector<BinnedED>>> scanComponentDists;
        std::map<std::string, std::map<std::string, std::vector<double>>> scanReactorRatios;
    };

    BuildResult BuildDatasets(
        FitInputs &in,
        const OutputDirs &dirs,
        const DatasetSetup &setup,
        const SystSetup &syst,
        const OscParams &osc,
        const std::unordered_map<int, double> &indexDistance,
        const DatasetParsInfo &dsParsInfo,
        const BuildDatasetsConfig &cfg = BuildDatasetsConfig());

    // Build one LLH per dataset, optionally overriding oscillated PDFs
    // with the scan-point-specific components from BuildResult::scanComponentDists.
    std::vector<BinnedNLLH> BuildScanPointLLHs(
        const BuildResult &datasets,
        const SystSetup &syst,
        const FitInputs &in,
        const PDFConfig &pdfConfig,
        size_t scanPointIndex);

} // namespace antinufit
