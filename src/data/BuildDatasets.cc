#include <BuildDatasets.hh>
#include <BuildComponent.hh>
#include <BuildDist.hh>
#include <LoadEvents.hh>
#include <RateReporting.hh>

#include <iostream>
#include <memory>
#include <stdexcept>

namespace antinufit
{
    void ApplyReactorRateScales(FitParameters &p,
                                const PreparedDatasetMap &datasets)
    {
        // LoadAnalysisConfig rejects reactor components shared across datasets
        for (PreparedDatasetMap::const_iterator dataset = datasets.begin();
             dataset != datasets.end(); ++dataset)
        {
            for (ParameterDict::const_iterator entry = dataset->second.reactorRatio.begin();
                 entry != dataset->second.reactorRatio.end(); ++entry)
            {
                const std::string &name = entry->first;
                const double scale = entry->second;
                p.noms.at(name) *= scale;
                p.mins.at(name) *= scale;
                p.maxs.at(name) *= scale;
                ParameterDict::iterator prior = p.constraintMeans.find(name);
                if (prior != p.constraintMeans.end())
                {
                    prior->second *= scale;
                    p.constraintSigmas.at(name) *= scale;
                }
            }
            for (ParameterDict::const_iterator entry = dataset->second.reactorFakeRatio.begin();
                 entry != dataset->second.reactorFakeRatio.end(); ++entry)
                p.fdValues.at(entry->first) *= entry->second;
        }
    }

    // Prepare likelihood inputs and output distributions for each dataset
    PreparedDatasetMap BuildDatasets(
        const FitConfig &fitConfig, const FitParameters &p,
        const DatasetSetup &setup,
        const SystsByDataset &syst,
        const SystConfig &systConfig,
        const OscParams &osc,
        const std::unordered_map<int, double> &indexDistance,
        const DatasetParameterMap &datasetPars)
    {
        PreparedDatasetMap out;

        for (DSMap::const_iterator dsIt = setup.dsPDFMap.begin(); dsIt != setup.dsPDFMap.end(); ++dsIt)
        {
            const std::string &dsName = dsIt->first;
            PreparedDataset &dataset = out[dsName];

            std::cout << "Preparing dataset: " << dsName << std::endl;
            std::vector<BinnedED> fakePdfs;
            if (fitConfig.GetFakeData())
                fakePdfs.reserve(dsIt->second.size());

            // Stage 1: load MC and build normalised component PDFs and reactor fractions
            for (EvMap::const_iterator evIt = dsIt->second.begin(); evIt != dsIt->second.end(); ++evIt)
            {
                const std::string &pdfName = evIt->first;
                const EventConfig &ev = evIt->second;
                std::cout << "Building distribution for " << pdfName << std::endl;

                std::unique_ptr<DataSet> dataSet = LoadEvents(ev.GetPrunedPath());
                if (!dataSet)
                    throw std::runtime_error("Cannot load MC for " + pdfName + " in " + dsName +
                                             ": " + ev.GetPrunedPath());

                dataset.pdfGroups.push_back(ev.GetGroup());
                const int numDimensions = ev.GetNumDimensions();
                BinnedED componentDist;
                BinnedED fakeComponentDist;
                double reactorRatio = 1.0;
                double reactorFakeRatio = 1.0;
                if (ev.GetOscillated())
                {
                    const double theta12SinSq = osc.theta12_sinsq();
                    componentDist = BuildDist::BuildOscillatedDist(
                        pdfName, numDimensions, setup.pdfConfig, dataSet.get(),
                        osc.deltam21, theta12SinSq, indexDistance, reactorRatio);

                    if (fitConfig.GetFakeData())
                    {
                        const double fakeTheta12SinSq = Theta12ToSinSq(
                            osc.theta12name, p.fdValues.at(osc.theta12name));
                        fakeComponentDist = BuildDist::BuildOscillatedDist(
                            pdfName, numDimensions, setup.pdfConfig, dataSet.get(),
                            p.fdValues.at("deltam21"), fakeTheta12SinSq,
                            indexDistance, reactorFakeRatio);
                    }

                }
                else if (ev.GetFlat())
                {
                    componentDist = BuildDist::BuildFlatDist(
                        pdfName, numDimensions, setup.pdfConfig);
                    if (fitConfig.GetFakeData())
                        fakeComponentDist = componentDist;
                }
                else
                {
                    componentDist = BuildDist::Build(
                        pdfName, numDimensions, setup.pdfConfig, dataSet.get());
                    if (fitConfig.GetFakeData())
                        fakeComponentDist = componentDist;
                }

                componentDist.AddPadding();
                if (fitConfig.GetFakeData())
                    fakeComponentDist.AddPadding();

                const int generated = static_cast<int>(componentDist.Integral());
                if (componentDist.Integral())
                    componentDist.Normalise();
                if (fitConfig.GetFakeData() && fakeComponentDist.Integral())
                    fakeComponentDist.Normalise();
                if (fitConfig.GetFakeData())
                    fakePdfs.push_back(fakeComponentDist);

                if (ev.GetOscillated())
                {
                    dataset.reactorRatio[pdfName] = reactorRatio;
                    if (fitConfig.GetFakeData())
                        dataset.reactorFakeRatio[pdfName] = reactorFakeRatio;
                }

                dataset.genRates.push_back(generated);
                dataset.pdfs.push_back(componentDist);
                dataset.normFittingStatuses.push_back(INDIRECT);

            }

            // Stage 2: use the PDFs and configured rates to build Asimov and fake data
            dataset.asimov = BinnedED("asimov", setup.systAxes);
            dataset.asimov.SetObservables(setup.dataObs);
            if (fitConfig.GetFakeData())
            {
                dataset.fakeDataset = BinnedED("fake_dataset", setup.systAxes);
                dataset.fakeDataset.SetObservables(setup.dataObs);
            }

            std::vector<double> nominalRateScales;
            std::vector<double> fakeRateScales;
            nominalRateScales.reserve(dataset.pdfs.size());
            if (fitConfig.GetFakeData())
                fakeRateScales.reserve(dataset.pdfs.size());
            for (const BinnedED &pdf : dataset.pdfs)
            {
                const std::map<std::string, double>::const_iterator ratio =
                    dataset.reactorRatio.find(pdf.GetName());
                nominalRateScales.push_back(ratio == dataset.reactorRatio.end() ? 1.0 : ratio->second);
                if (fitConfig.GetFakeData())
                {
                    const std::map<std::string, double>::const_iterator fakeRatio =
                        dataset.reactorFakeRatio.find(pdf.GetName());
                    fakeRateScales.push_back(fakeRatio == dataset.reactorFakeRatio.end()
                                                 ? 1.0 : fakeRatio->second);
                }
            }
            std::vector<BinnedED> nominalComponents = BuildComponents(
                dataset.pdfs, dataset.pdfGroups, syst, systConfig, dsName,
                p.noms, nominalRateScales, dataset.asimov, setup.dataObs);
            dataset.nominalNonBufferRates = GetComponentRates(
                dataset.pdfs, nominalComponents, setup.pdfConfig, setup.dataObs);
            if (fitConfig.GetSaveOutputs())
                dataset.nominalComponents = std::move(nominalComponents);

            if (fitConfig.GetFakeData())
            {
                std::vector<BinnedED> fakeComponents = BuildComponents(
                    fakePdfs, dataset.pdfGroups, syst, systConfig, dsName,
                    p.fdValues, fakeRateScales, dataset.fakeDataset, setup.dataObs);
                if (fitConfig.GetSaveOutputs())
                    dataset.fakeComponents = std::move(fakeComponents);
            }

            // Stage 3: select the distribution that the likelihood will treat as data
            if (fitConfig.GetAsimov())
                dataset.dataDist = dataset.asimov;
            else if (fitConfig.GetFakeData())
                dataset.dataDist = dataset.fakeDataset;
            else
                dataset.dataDist = LoadDataDist(setup.dataPath.at(dsName), setup.pdfConfig);

            // Stage 4: initialise this dataset's parameters with reactor rate scaling
            dataset.initialValues = InitialiseDatasetParams(p.noms, datasetPars, dsName);
            for (EvMap::const_iterator component = dsIt->second.begin();
                 component != dsIt->second.end(); ++component)
                if (component->second.GetOscillated())
                    dataset.initialValues.at(component->first) *= dataset.reactorRatio.at(component->first);
            std::cout << "Prepared dataset: " << dsName << std::endl << std::endl;
        }

        return out;
    }

} // namespace antinufit
