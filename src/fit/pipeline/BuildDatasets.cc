#include <BuildDatasets.hh>
#include <ApplySystematics.hh>
#include <BuildComponent.hh>
#include <Projection.hh>
#include <LoadDataset.hh>
#include <RateReporting.hh>

#include <IO.h>
#include <Shape.h>
#include <Event.h>
#include <algorithm>
#include <cmath>

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
            SystSetup &syst,
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

            std::string oscillationName;
            bool floatingOscillations = false;
            const std::vector<std::string> groups = ev.GetGroup();
            for (const auto &entry : syst.type)
            {
                if (entry.second != "FusedOsc" && entry.second.rfind("FusedOsc:", 0) != 0)
                    continue;
                floatingOscillations = true;
                const auto &datasets = syst.dataSets.at(entry.first);
                if (std::find(datasets.begin(), datasets.end(), dsName) == datasets.end())
                    continue;
                const std::string &group = syst.group.at(entry.first);
                if (group.empty())
                    throw std::runtime_error("FusedOsc needs an explicit group: " + entry.first);
                if (std::find(groups.begin(), groups.end(), group) == groups.end())
                    continue;
                if (!ev.GetOscillated() || groups.front() != group || ev.GetNumDimensions() != 3 ||
                    !oscillationName.empty() || syst.byDataset[dsName].count(entry.first))
                    throw std::runtime_error("FusedOsc must be first and belong to one oscillated 3D component: " + entry.first);
                oscillationName = entry.first;
            }
            if (floatingOscillations && ev.GetOscillated() && oscillationName.empty())
                throw std::runtime_error("No FusedOsc group for oscillated component " + pdfName);

            // Get MC events from file
            std::unique_ptr<DataSet> dataSet = LoadMCDataset(ev.GetPrunedPath());
            if (!dataSet)
            {
                if (!oscillationName.empty())
                    throw std::runtime_error("Cannot load reactor MC for " + pdfName);
                return;
            }
            model.pdfGroups.push_back(ev.GetGroup());

            const int numDimensions = ev.GetNumDimensions();
            BinnedED reactorMC;
            if (!oscillationName.empty())
            {
                const auto &params = syst.paramNames.at(oscillationName);
                const std::string &type = syst.type.at(oscillationName);
                const auto &obs = syst.distObs.at(oscillationName);
                if (params.size() != 2 || params.front() != "deltam21" ||
                    obs != std::vector<std::string>{"energy"} || syst.transObs.at(oscillationName) != obs ||
                    type != "FusedOsc")
                    throw std::runtime_error("Invalid FusedOsc parameters, observables or mode: " + oscillationName);
                for (const auto &entry : syst.byDataset[dsName])
                    if (syst.group.at(entry.first).empty() || syst.group.at(entry.first) == syst.group.at(oscillationName))
                        throw std::runtime_error("FusedOsc must have its own group: " + oscillationName);

                // Check reactor identities before histogram filling can clamp them into another bin.
                const BinAxis reactorAxis = DistBuilder::BuildAxes(setup.pdfConfig, 3).GetAxis(1);
                for (size_t i = 0; i < dataSet->GetNEntries(); ++i)
                {
                    const double reactor = dataSet->GetEntry(i).GetDatum("reactorIndex");
                    if (!std::isfinite(reactor) || reactor < reactorAxis.GetMin() || reactor >= reactorAxis.GetMax() ||
                        std::fabs(reactor - std::round(reactor)) > 1e-6 || !indexDistance.count(int(std::lround(reactor))))
                        throw std::runtime_error("Invalid reactor index or missing baseline in " + pdfName + " event " + std::to_string(i));
                }
                reactorMC = DistBuilder::Build(pdfName, 3, setup.pdfConfig, dataSet.get());
                std::unique_ptr<Shape> shape(new Shape(oscillationName));
                shape->SetShapeFunction(MakeReactorOscillationShape(reactorMC, indexDistance, params[1]), params);
                for (const std::string &param : params)
                    shape->SetParameter(param, in.p.noms.at(param));
                shape->SetAxes(DistBuilder::BuildAxes(setup.pdfConfig, 1));
                shape->SetDistributionObs(obs);
                shape->SetTransformationObs(syst.transObs.at(oscillationName));
                syst.byDataset[dsName][oscillationName] = shape.get();
                syst.owned.push_back(std::move(shape));
            }

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
                cfg.component, oscillationName.empty() ? nullptr : &reactorMC);

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

            if (oscillationName.empty() || in.run.isFakeData)
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
            out.nominalNonBufferRates[dsName][pdfName] += IntegrateBufferedBins(
                asimovComp, setup.pdfConfig, setup.dataObs);

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
                                const PDFConfig &pdfConfig,
                                const SystSetup &syst,
                                const std::string &dsName,
                                DatasetModel &model,
                                const BinnedED &dataDist)
        {
            BinnedNLLH &lh = llhs.emplace_back();
            for (std::vector<std::string>::const_iterator axisIt = pdfConfig.GetAxisNames().begin();
                 axisIt != pdfConfig.GetAxisNames().end(); ++axisIt)
            {
                if (!pdfConfig.HasLLHBufferBins(*axisIt))
                    continue;

                const std::pair<int, int> bufferBins = pdfConfig.GetLLHBufferBins(*axisIt);
                lh.SetBuffer(*axisIt, bufferBins.first, bufferBins.second);
            }

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
        SystSetup &syst,
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

            for (const auto &entry : syst.type)
                if ((entry.second == "FusedOsc" || entry.second.rfind("FusedOsc:", 0) == 0) &&
                    std::find(syst.dataSets.at(entry.first).begin(), syst.dataSets.at(entry.first).end(), dsName) !=
                        syst.dataSets.at(entry.first).end() && !syst.byDataset[dsName].count(entry.first))
                    throw std::runtime_error("Unused FusedOsc systematic " + entry.first + " in " + dsName);

            std::cout << std::endl;
            SaveFinalDatasetDists(in, dirs, dsName, model);
            std::cout << std::endl;

            const BinnedED dataDist = ChooseDataDist(in, setup, dsName, model);
            out.dataDists[dsName] = dataDist;

            BinnedNLLH &lh = BuildDatasetLLH(out.llhs, in, setup.pdfConfig, syst, dsName, model, dataDist);

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
        const PDFConfig &pdfConfig,
        size_t scanPointIndex)
    {
        std::vector<BinnedNLLH> llhs;
        llhs.reserve(datasets.models.size());

        for (DatasetModelMap::const_iterator dsIt = datasets.models.begin(); dsIt != datasets.models.end(); ++dsIt)
        {
            const std::string &dsName = dsIt->first;
            const DatasetModel &model = dsIt->second;

            BinnedNLLH &lh = llhs.emplace_back();
            for (std::vector<std::string>::const_iterator axisIt = pdfConfig.GetAxisNames().begin();
                 axisIt != pdfConfig.GetAxisNames().end(); ++axisIt)
            {
                if (!pdfConfig.HasLLHBufferBins(*axisIt))
                    continue;

                const std::pair<int, int> bufferBins = pdfConfig.GetLLHBufferBins(*axisIt);
                lh.SetBuffer(*axisIt, bufferBins.first, bufferBins.second);
            }
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
