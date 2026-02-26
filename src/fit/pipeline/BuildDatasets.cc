#include <BuildDatasets.hh>
#include <ApplySystematics.hh>
#include <BuildComponent.hh>
#include <Projection.hh>
#include <LoadDataset.hh>

#include <IO.h>

namespace antinufit
{
    namespace
    {
        void InitialiseDatasetModel(DatasetModel &model, const DatasetSetup &setup)
        {
            // Create the empty full dist
            model.asimov = BinnedED("asimov", setup.systAxes);
            model.asimov.SetObservables(setup.dataObs);

            // And an empty fake data dist
            model.fakeDataset = BinnedED("fake_dataset", setup.systAxes);
            model.fakeDataset.SetObservables(setup.dataObs);
        }

        void BuildAndAccumulateComponent(
            FitInputs &in,
            const OutputDirs &dirs,
            const DatasetSetup &setup,
            const SystSetup &syst,
            const BuildDatasetsConfig &cfg,
            const OscParams &osc,
            const std::unordered_map<int, double> &indexDistance,
            const std::string &dsName,
            const std::string &pdfName,
            const EventConfig &ev,
            DatasetModel &model,
            std::map<std::string, double> &reactorRatio,
            std::map<std::string, double> &reactorRatioFD,
            BuildResult &out)
        {
            std::cout << "Building distribution for " << pdfName << std::endl;
            model.pdfGroups.push_back(ev.GetGroup());

            // Get MC events from file
            std::unique_ptr<DataSet> dataSet = LoadMCDataset(ev.GetPrunedPath());
            if (!dataSet)
                return;

            const int numDimensions = ev.GetNumDimensions();

            // Build the component dists (no systematics yet)
            ComponentBuildResult comp = BuildComponentDists(
                in,
                setup,
                ev,
                pdfName,
                numDimensions,
                dataSet.get(),
                osc,
                indexDistance,
                reactorRatio,
                reactorRatioFD,
                cfg.component);

            // Beeston-Barlow generated rate (pre-normalise count)
            model.genRates.push_back(comp.generated);

            // Store the normalised base pdf for LLH building
            model.pdfs.push_back(comp.dist);
            model.normFittingStatuses.push_back(INDIRECT);

            // Snapshot of the unscaled pdf (normalised, pre-syst, pre-rate scaling)
            BinnedED unscaledPDF = comp.unscaledPdf;

            // Working copies for asimov/fake building
            BinnedED dist = comp.dist;
            BinnedED fakeDataDist = comp.fakeDist;

            // Apply nominal systematics
            dist = ApplySystematicsByGroups(
                std::move(dist),
                model.pdfGroups.back(),
                syst.byDataset.at(dsName),
                syst.group,
                in.p.noms);

            fakeDataDist = ApplySystematicsByGroups(
                std::move(fakeDataDist),
                model.pdfGroups.back(),
                syst.byDataset.at(dsName),
                syst.group,
                in.p.fdValues);

            // Now scale the Asimov component by expected count, and also save pdf as a histo
            dist.Scale(in.p.noms.at(pdfName));
            fakeDataDist.Scale(in.p.fdValues.at(pdfName));

            const std::string baseFile = pdfName + "_" + dsName + ".root";

            BinnedED asimovComp = ProjectToDataObs(dist, model.asimov, setup.dataObs);
            model.asimov.Add(asimovComp);

            if (in.run.saveOutputs)
            {
                IO::SaveHistogram(asimovComp.GetHistogram(),
                                  dirs.asimovDistDir + "/" + baseFile,
                                  dist.GetName());

                BinnedED pdfToSave = ProjectToDataObs(unscaledPDF, model.asimov, setup.dataObs);
                IO::SaveHistogram(pdfToSave.GetHistogram(),
                                  dirs.pdfDir + "/" + baseFile,
                                  unscaledPDF.GetName());
            }

            if (in.run.isFakeData)
            {
                BinnedED fakeComp = ProjectToDataObs(fakeDataDist, model.asimov, setup.dataObs);
                model.fakeDataset.Add(fakeComp);

                if (in.run.saveOutputs)
                {
                    IO::SaveHistogram(fakeComp.GetHistogram(),
                                      dirs.fakedataDistDir + "/" + baseFile,
                                      fakeDataDist.GetName());
                }
            }

            if (!comp.scanDists.empty())
            {
                out.scanComponentDists[dsName][pdfName] = comp.scanDists;
                out.scanReactorRatios[dsName][pdfName] = comp.scanReactorRatios;
            }
        }

        void SaveFinalDatasetDists(const FitInputs &in,
                                   const OutputDirs &dirs,
                                   const std::string &dsName,
                                   const DatasetModel &model)
        {
            // Save combined histogram (final asimov dataset)
            if (in.run.saveOutputs)
            {
                IO::SaveHistogram(model.asimov.GetHistogram(), dirs.outDir + "/asimov_" + dsName + ".root", "asimov");
                if (in.run.isFakeData)
                {
                    IO::SaveHistogram(model.fakeDataset.GetHistogram(), dirs.outDir + "/fakedata_" + dsName + ".root", "fakedata");
                }
            }
        }

        BinnedED ChooseDataDist(const FitInputs &in,
                                const DatasetSetup &setup,
                                const std::string &dsName,
                                const DatasetModel &model)
        {
            if (in.run.isAsimov)
                return model.asimov;

            if (in.run.isFakeData)
                return model.fakeDataset;

            return LoadDataDist(setup.dataPath.at(dsName), setup.pdfConfig);
        }

        BinnedNLLH &BuildDatasetLLH(std::vector<BinnedNLLH> &llhs,
                                    const FitInputs &in,
                                    const SystSetup &syst,
                                    const std::string &dsName,
                                    DatasetModel &model,
                                    const BinnedED &dataDist)
        {
            BinnedNLLH &lh = llhs.emplace_back();
            lh.SetBuffer("energy", 8, 20);

            // Add our data
            lh.SetDataDist(dataDist);

            // Set whether or not to use Beeston Barlow
            lh.SetBarlowBeeston(in.run.beestonBarlowFlag);

            // Add the systematics and any prior constraints
            for (std::map<std::string, Systematic *>::const_iterator systIt = syst.byDataset.at(dsName).begin();
                 systIt != syst.byDataset.at(dsName).end(); ++systIt)
            {
                lh.AddSystematic(systIt->second, syst.group.at(systIt->first));
            }

            // Add our pdfs
            lh.AddPdfs(model.pdfs, model.pdfGroups, model.genRates, &model.normFittingStatuses);

            // Bring it all together
            lh.RegisterFitComponents();
            return lh;
        }
    } // namespace

    // Build the PDFs (scaled and distorted) and LLHs for one dataset
    BuildResult BuildDatasets(
        FitInputs &in,
        const OutputDirs &dirs,
        const DatasetSetup &setup,
        const SystSetup &syst,
        const OscParams &osc,
        const std::unordered_map<int, double> &indexDistance,
        const DatasetParsInfo &dsParsInfo,
        const BuildDatasetsConfig &cfg)
    {
        BuildResult out;
        out.llhs.reserve(setup.dsPDFMap.size());

        // Loop over datasets and build their distributions + LLHs
        for (DSMap::const_iterator dsIt = setup.dsPDFMap.begin(); dsIt != setup.dsPDFMap.end(); ++dsIt)
        {
            const std::string &dsName = dsIt->first;
            DatasetModel &model = out.models[dsName];

            std::cout << "Building Asimov for dataset: " << dsName << std::endl;
            InitialiseDatasetModel(model, setup);

            // Build each PDF component and accumulate into dataset products
            for (EvMap::const_iterator evIt = dsIt->second.begin(); evIt != dsIt->second.end(); ++evIt)
            {
                BuildAndAccumulateComponent(
                    in,
                    dirs,
                    setup,
                    syst,
                    cfg,
                    osc,
                    indexDistance,
                    dsName,
                    evIt->first,
                    evIt->second,
                    model,
                    out.reactorRatio,
                    out.reactorRatioFD,
                    out);
            }

            std::cout << std::endl;
            SaveFinalDatasetDists(in, dirs, dsName, model);
            std::cout << std::endl;

            const BinnedED dataDist = ChooseDataDist(in, setup, dsName, model);
            out.dataDists[dsName] = dataDist;

            BinnedNLLH &lh = BuildDatasetLLH(out.llhs, in, syst, dsName, model, dataDist);

            out.initialByDataset[dsName] = InitialiseDatasetParams(in.p.noms, dsParsInfo.datasetPars, dsName);

            // Set to these initial values
            lh.SetParameters(out.initialByDataset.at(dsName));

            std::cout << "Made LLH for Dataset: " << dsName << std::endl
                      << std::endl;
        }

        return out;
    }

    std::vector<BinnedNLLH> BuildScanPointLLHs(
        const BuildResult &datasets,
        const SystSetup &syst,
        const FitInputs &in,
        size_t scanPointIndex)
    {
        std::vector<BinnedNLLH> llhs;
        llhs.reserve(datasets.models.size());

        for (DatasetModelMap::const_iterator dsIt = datasets.models.begin(); dsIt != datasets.models.end(); ++dsIt)
        {
            const std::string &dsName = dsIt->first;
            const DatasetModel &model = dsIt->second;

            BinnedNLLH &lh = llhs.emplace_back();
            lh.SetBuffer("energy", 8, 20);
            lh.SetDataDist(datasets.dataDists.at(dsName));
            lh.SetBarlowBeeston(in.run.beestonBarlowFlag);

            for (std::map<std::string, Systematic *>::const_iterator systIt = syst.byDataset.at(dsName).begin();
                 systIt != syst.byDataset.at(dsName).end(); ++systIt)
            {
                lh.AddSystematic(systIt->second, syst.group.at(systIt->first));
            }

            std::vector<BinnedED> pdfs = model.pdfs;
            std::map<std::string, std::map<std::string, std::vector<BinnedED>>>::const_iterator dsScanIt =
                datasets.scanComponentDists.find(dsName);

            if (dsScanIt != datasets.scanComponentDists.end())
            {
                for (std::map<std::string, std::vector<BinnedED>>::const_iterator pdfIt = dsScanIt->second.begin();
                     pdfIt != dsScanIt->second.end(); ++pdfIt)
                {
                    if (scanPointIndex >= pdfIt->second.size())
                        continue;

                    for (size_t iPdf = 0; iPdf < pdfs.size(); ++iPdf)
                    {
                        if (pdfs[iPdf].GetName() == pdfIt->first)
                        {
                            pdfs[iPdf] = pdfIt->second[scanPointIndex];
                            break;
                        }
                    }
                }
            }

            lh.AddPdfs(pdfs, model.pdfGroups, model.genRates, &model.normFittingStatuses);
            lh.RegisterFitComponents();
            lh.SetParameters(datasets.initialByDataset.at(dsName));
        }

        return llhs;
    }

} // namespace antinufit
