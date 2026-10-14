// Antinu headers
#include <PDFConfigLoader.hh>
#include <DistBuilder.hh>
#include <DistTools.h>
#include <FitConfigLoader.hh>
#include <EventConfigLoader.hh>
#include <SystConfigLoader.hh>
#include <SystFactory.hh>
#include <OscGridConfigLoader.hh>
#include <OscWeight1D.hh>
#include <Utilities.hh>

// OXO headers
#include <ROOTNtuple.h>
#include <BinnedNLLH.h>
#include <StatisticSum.h>
#include <IO.h>
#include <Rand.h>
#include <MetropolisSampler.h>
#include <MCMC.h>

// ROOT headers
#include <TH1D.h>
#include <TKey.h>
#include <TStopwatch.h>

// c++ headers
#include <sys/stat.h>

using namespace antinufit;

namespace
{
  bool IsFusedOsc(const std::string &type_)
  {
    return type_.rfind("FusedOsc", 0) == 0;
  }

  // Apply every systematic of each group in groups_, group by group in that order, as BinnedNLLH does.
  // Systematic parameters are first set from params_. The integral change is kept (&norm), so the result
  // carries oscillation and edge losses exactly as the fit's working PDF does.
  BinnedED ApplyGroups(BinnedED dist_, const std::vector<std::string> &groups_,
                       std::map<std::string, Systematic *> &systs_,
                       std::map<std::string, std::string> &systGroup_,
                       const ParameterDict &params_, std::vector<double> *integrals_ = nullptr)
  {
    for (size_t iGroup = 0; iGroup < groups_.size(); iGroup++)
    {
      for (std::map<std::string, Systematic *>::iterator systIt = systs_.begin(); systIt != systs_.end(); ++systIt)
      {
        if (systGroup_[systIt->first] != groups_.at(iGroup))
          continue;
        std::set<std::string> systParamNames = systIt->second->GetParameterNames();
        for (std::set<std::string>::iterator parIt = systParamNames.begin(); parIt != systParamNames.end(); ++parIt)
          systIt->second->SetParameter(*parIt, params_.at(*parIt));
        systIt->second->Construct();
        double norm; // Passing &norm stops OXSX re-normalising the result; the value written into it is not used
        dist_ = systIt->second->operator()(dist_, &norm);
      }
      if (integrals_)
        integrals_->push_back(dist_.Integral());
    }
    return dist_;
  }
}

void mcmc(const std::string &fitConfigFile_,
          const std::string &evConfigFile_,
          const std::string &pdfConfigFile_,
          const std::string &systConfigFile_,
          const std::string &oscGridConfigFile_)
{
  Rand::SetSeed(0);

  // Load up the fit configuration information
  FitConfig fitConfig;
  FitConfigLoader fitLoader(fitConfigFile_);
  fitConfig = fitLoader.LoadActive();
  bool isAsimov = fitConfig.GetAsimov();
  bool isFakeData = fitConfig.GetFakeData();
  bool beestonBarlowFlag = fitConfig.GetBeestonBarlow();
  int nsteps = fitConfig.GetIterations();
  int burnin = fitConfig.GetBurnIn();
  std::string outDir = fitConfig.GetOutDir();
  ParameterDict constrMeans = fitConfig.GetConstrMeans();
  ParameterDict constrSigmas = fitConfig.GetConstrSigmas();
  ParameterDict mins = fitConfig.GetMinima();
  ParameterDict maxs = fitConfig.GetMaxima();
  ParameterDict noms = fitConfig.GetNominals();
  ParameterDict sigmas = fitConfig.GetSigmas();
  ParameterDict nbins = fitConfig.GetNBins();
  ParameterDict constrRatioMeans = fitConfig.GetConstrRatioMeans();
  ParameterDict constrRatioSigmas = fitConfig.GetConstrRatioSigmas();
  std::map<std::string, std::string> constrRatioParName = fitConfig.GetConstrRatioParName();
  ParameterDict constrCorrs = fitConfig.GetConstrCorrs();
  std::map<std::string, std::string> constrCorrParName = fitConfig.GetConstrCorrParName();
  ParameterDict fdValues = fitConfig.GetFakeDataVals();
  std::map<std::string, bool> fixedPars = fitConfig.GetFixPars();
  double sigmaScale = fitConfig.GetSigmaScale();

  std::string pdfDir = outDir + "/unscaled_pdfs";
  std::string asimovDistDir = outDir + "/asimov_dists";
  std::string fakedataDistDir = outDir + "/fakedata_dists";
  std::string postfitDistDir = outDir + "/postfit_dists";
  std::string projDir1D = outDir + "/1dlhproj";
  std::string projDir2D = outDir + "/2dlhproj";

  // Create output directories
  struct stat st = {0};
  for (const std::string &dir : {outDir, pdfDir, asimovDistDir, fakedataDistDir, postfitDistDir, projDir1D, projDir2D})
    if (stat(dir.c_str(), &st) == -1)
      mkdir(dir.c_str(), 0700);

  // Load up all the event types we want to contribute
  typedef std::map<std::string, EventConfig> EvMap;
  typedef std::map<std::string, std::map<std::string, EventConfig>> DSMap;
  EventConfigLoader evLoader(evConfigFile_);
  DSMap dsPDFMap = evLoader.LoadActive();
  std::map<std::string, std::string> dataPath = evLoader.GetPrunedDataPaths();
  std::map<std::string, BinnedED> dataDists;

  // Load up the PDF information (skeleton axis details, rather than the distributions themselves)
  PDFConfigLoader pdfLoader(pdfConfigFile_);
  PDFConfig pdfConfig = pdfLoader.Load();
  std::vector<std::string> dataObs = pdfConfig.GetDataBranchNames();
  ObsSet dataObsSet(dataObs);
  AxisCollection systAxes = DistBuilder::BuildAxes(pdfConfig, dataObs.size());

  // Load up the systematics
  SystConfigLoader systLoader(systConfigFile_);
  SystConfig systConfig = systLoader.LoadActive();
  std::map<std::string, std::vector<std::string>> systParamNames = systConfig.GetParamNames();
  std::map<std::string, std::string> systGroup = systConfig.GetGroup();
  std::map<std::string, std::string> systType = systConfig.GetType();
  std::map<std::string, std::vector<std::string>> systDistObs = systConfig.GetDistObs();
  std::map<std::string, std::vector<std::string>> systTransObs = systConfig.GetTransObs();
  std::map<std::string, std::vector<std::string>> systDataSets = systConfig.GetDataSets();
  std::map<std::string, std::vector<std::string>> systParDataSets = systConfig.GetParDataSets();
  std::vector<std::string> fullParamNameVec;

  // No systematic may sit in the global group: OXSX applies it before all listed groups, which would put
  // it ahead of the oscillation
  for (std::map<std::string, std::string>::iterator systIt = systGroup.begin(); systIt != systGroup.end(); ++systIt)
    if (systIt->second == "")
      throw std::runtime_error("ERROR: Systematic " + systIt->first + " has no group. Every systematic must name a group.");

  // Reactor baselines, and the oscillation grids only if a systematic reads them
  OscGridConfigLoader oscGridLoader(oscGridConfigFile_);
  OscGridConfig oscGridConfig = oscGridLoader.Load();
  std::unordered_map<int, double> indexDistance = LoadIndexDistanceMap(oscGridConfig.GetReactorsJsonFile());
  std::map<int, OscGrid *> oscGridMap;
  for (std::map<std::string, std::string>::iterator systIt = systType.begin(); systIt != systType.end(); ++systIt)
  {
    if (systIt->second != "Shape:OscProbGrid" || oscGridMap.size())
      continue;
    for (std::unordered_map<int, double>::iterator it = indexDistance.begin(); it != indexDistance.end(); ++it)
      oscGridMap[it->first] = new OscGrid(oscGridConfig.GetFilename() + "_" + std::to_string(it->first) + ".root");
  }

  // Loop over systematics and declare each one. A FusedOsc systematic holds one reactor PDF's MC, so it
  // gets its own instance per dataset; every other systematic is shared between datasets.
  std::map<std::string, std::map<std::string, Systematic *>> systMap;
  std::map<std::string, std::map<std::string, OscWeight1D *>> fusedOscMap;
  for (std::map<std::string, std::string>::iterator systIt = systType.begin(); systIt != systType.end(); ++systIt)
  {
    // Check any parameter defined for this systematic was declared in the fit config
    for (size_t iParam = 0; iParam < systParamNames[systIt->first].size(); iParam++)
    {
      if (noms.find(systParamNames[systIt->first].at(iParam)) == noms.end())
        throw std::runtime_error("ERROR: Syst config defines systematic " + systIt->first + " with parameter " +
                                 systParamNames[systIt->first].at(iParam) + ", but this parameter is not defined in the fit config");
    }
    fullParamNameVec.insert(fullParamNameVec.end(), systParamNames[systIt->first].begin(), systParamNames[systIt->first].end());

    for (std::vector<std::string>::iterator dsIt = systDataSets[systIt->first].begin(); dsIt != systDataSets[systIt->first].end(); ++dsIt)
      if (dsPDFMap.find(*dsIt) == dsPDFMap.end())
        throw std::runtime_error("ERROR: SystConfig defines systematic " + systIt->first + " applying to dataset " + *dsIt + " that is not defined in EventConfig");

    AxisCollection axes = DistBuilder::BuildAxes(pdfConfig, systDistObs[systIt->first].size());
    Systematic *shared = nullptr;
    for (std::vector<std::string>::iterator dsIt = systDataSets[systIt->first].begin(); dsIt != systDataSets[systIt->first].end(); ++dsIt)
    {
      Systematic *syst = shared;
      if (!syst)
      {
        syst = SystFactory::New(systIt->first, systType[systIt->first], systParamNames[systIt->first], noms, oscGridMap, indexDistance);
        syst->SetAxes(axes);
        // The "dimensions" the systematic applies to
        syst->SetTransformationObs(systTransObs[systIt->first]);
        // All the "dimensions" of the dataset
        syst->SetDistributionObs(systDistObs[systIt->first]);
        // A FusedOsc systematic is constructed once its reactor MC is attached
        if (!IsFusedOsc(systType[systIt->first]))
        {
          syst->Construct();
          shared = syst;
        }
      }
      systMap[*dsIt][systIt->first] = syst;
      if (IsFusedOsc(systType[systIt->first]))
        fusedOscMap[*dsIt][systIt->first] = dynamic_cast<OscWeight1D *>(syst);
    }
  }

  // Every oscillated event type must list a FusedOsc group first, and that group must hold only FusedOsc
  // systematics. diag(P_ee) and the energy response do not commute, so nothing may act on the reactor PDF
  // before the oscillation. No other event type may list that group.
  std::map<std::string, std::map<std::string, OscWeight1D *>> reactorOsc; // [dataset][event type]
  for (DSMap::iterator dsIt = dsPDFMap.begin(); dsIt != dsPDFMap.end(); ++dsIt)
  {
    std::map<std::string, OscWeight1D *> fusedByGroup;
    for (std::map<std::string, OscWeight1D *>::iterator oscIt = fusedOscMap[dsIt->first].begin(); oscIt != fusedOscMap[dsIt->first].end(); ++oscIt)
    {
      if (fusedByGroup.count(systGroup[oscIt->first]))
        throw std::runtime_error("ERROR: Group " + systGroup[oscIt->first] + " holds more than one FusedOsc systematic in dataset " + dsIt->first);
      fusedByGroup[systGroup[oscIt->first]] = oscIt->second;
    }
    for (std::map<std::string, Systematic *>::iterator systIt = systMap[dsIt->first].begin(); systIt != systMap[dsIt->first].end(); ++systIt)
      if (fusedByGroup.count(systGroup[systIt->first]) && !IsFusedOsc(systType[systIt->first]))
        throw std::runtime_error("ERROR: Systematic " + systIt->first + " shares group " + systGroup[systIt->first] + " with a FusedOsc systematic; that group may hold only the oscillation");

    std::set<OscWeight1D *> used;
    for (EvMap::iterator evIt = dsIt->second.begin(); evIt != dsIt->second.end(); ++evIt)
    {
      std::vector<std::string> groups = evIt->second.GetGroup();
      for (size_t iGroup = 0; iGroup < groups.size(); iGroup++)
      {
        if (!fusedByGroup.count(groups.at(iGroup)))
          continue;
        if (!evIt->second.GetOscillated() || iGroup != 0)
          throw std::runtime_error("ERROR: Event type " + evIt->first + " in dataset " + dsIt->first + " lists FusedOsc group " + groups.at(iGroup) +
                                   "; only an oscillated event type may, and only as its first group");
      }
      if (!evIt->second.GetOscillated())
        continue;
      if (groups.empty() || !fusedByGroup.count(groups.front()))
        throw std::runtime_error("ERROR: Oscillated event type " + evIt->first + " in dataset " + dsIt->first + " must list a FusedOsc group first in groups");
      if (evIt->second.GetNumDimensions() != 3)
        throw std::runtime_error("ERROR: Oscillated event type " + evIt->first + " must have dimensions = 3 (energy, reactor, nu_energy)");
      OscWeight1D *osc = fusedByGroup[groups.front()];
      if (used.count(osc))
        throw std::runtime_error("ERROR: Two oscillated event types in dataset " + dsIt->first + " share FusedOsc group " + groups.front());
      used.insert(osc);
      reactorOsc[dsIt->first][evIt->first] = osc;
    }
    for (std::map<std::string, OscWeight1D *>::iterator oscIt = fusedOscMap[dsIt->first].begin(); oscIt != fusedOscMap[dsIt->first].end(); ++oscIt)
      if (!used.count(oscIt->second))
        throw std::runtime_error("ERROR: FusedOsc systematic " + oscIt->first + " is not used by any oscillated event type in dataset " + dsIt->first);
  }

  // A parameter could have been defined in the fit config but isn't associated with a pdf or systematic
  // If so ignore that parameter
  ParameterDict::iterator parIt = noms.begin();
  size_t numParams = noms.size();
  std::map<std::string, std::vector<std::string>> datasetPars;
  for (size_t iParam = 0; iParam < numParams; iParam++)
  {
    // We will check if the parameter is a pdf or systematic for any dataset.
    // We'll also use 'datasetPars' to keep track of which datasets each parameter applies to
    bool isPDF = false;
    bool isSyst = false;
    bool iterate = true;
    for (DSMap::iterator dsIt = dsPDFMap.begin(); dsIt != dsPDFMap.end(); ++dsIt)
    {
      if (dsIt->second.find(parIt->first) != dsIt->second.end())
      {
        isPDF = true;
        datasetPars[parIt->first].push_back(dsIt->first);
      }
    }
    if (!isPDF && systParDataSets.find(parIt->first) != systParDataSets.end())
    {
      isSyst = true;
      datasetPars[parIt->first] = systParDataSets[parIt->first];
    }

    if (!isPDF && !isSyst)
    {
      std::cout << parIt->first << " parameter defined in fit config but not in syst or event config. It will be ignored." << std::endl;
      constrSigmas.erase(parIt->first);
      constrMeans.erase(parIt->first);
      constrRatioMeans.erase(parIt->first);
      constrRatioSigmas.erase(parIt->first);
      constrRatioParName.erase(parIt->first);
      constrCorrs.erase(parIt->first);
      constrCorrParName.erase(parIt->first);
      mins.erase(parIt->first);
      maxs.erase(parIt->first);
      sigmas.erase(parIt->first);
      nbins.erase(parIt->first);
      fixedPars.erase(parIt->first);
      parIt = noms.erase(parIt);
      if (parIt != noms.begin())
        parIt--;
      else
        iterate = false;
    }
    if (iterate)
      parIt++;
  }
  std::cout << std::endl;

  PrintParams(mins, maxs, noms, constrMeans, constrSigmas, constrRatioMeans, constrRatioSigmas, constrRatioParName, constrCorrs, constrCorrParName, datasetPars, fixedPars);

  // Create the individual PDFs and Asimov components, for each dataset, and make the component LLH objects
  std::map<std::string, std::vector<BinnedED>> pdfMap;
  std::map<std::string, std::vector<std::vector<std::string>>> pdfGroups;
  std::vector<BinnedNLLH> testStats;
  testStats.reserve(dsPDFMap.size());
  std::map<std::string, ParameterDict> parameterValues;
  std::map<std::string, std::vector<int>> genRates;
  std::map<std::string, std::vector<NormFittingStatus> *> normFittingStatuses;

  for (DSMap::iterator dsIt = dsPDFMap.begin(); dsIt != dsPDFMap.end(); ++dsIt)
  {
    std::cout << "Building Asimov for dataset: " << dsIt->first << std::endl;
    normFittingStatuses[dsIt->first] = new std::vector<NormFittingStatus>;

    // Create the empty full dist
    BinnedED asimov = BinnedED("asimov", systAxes);
    asimov.SetObservables(dataObs);
    // And an empty fake data dist
    BinnedED fakeDataset = BinnedED("fake_dataset", systAxes);
    fakeDataset.SetObservables(dataObs);

    // Build each of the PDFs, scale them to the correct size
    for (EvMap::iterator evIt = dsIt->second.begin(); evIt != dsIt->second.end(); ++evIt)
    {
      std::cout << "Building distribution for " << evIt->first << std::endl;

      // Get MC events from file
      DataSet *dataSet;
      try
      {
        std::cout << "Loading " << evIt->second.GetPrunedPath() << std::endl;
        dataSet = new ROOTNtuple(evIt->second.GetPrunedPath(), "pruned");
      }
      catch (const IOError &e_)
      {
        if (evIt->second.GetOscillated())
          throw;
        std::cout << "Warning: skipping " << evIt->first << " couldn't open files:\n\t" << e_.what() << std::endl;
        continue;
      }
      pdfGroups[dsIt->first].push_back(evIt->second.GetGroup());

      // Build distribution of those events
      BinnedED dist;
      int num_dimensions = evIt->second.GetNumDimensions();
      if (evIt->second.GetOscillated())
      {
        // The reactor PDF is registered in 1D; the 3D (energy, reactor, nu_energy) MC is handed to its
        // FusedOsc systematic as data
        OscWeight1D *osc = reactorOsc[dsIt->first][evIt->first];
        AxisCollection axes3d = DistBuilder::BuildAxes(pdfConfig, 3);
        size_t nClamped = OscWeight1D::CheckEvents(dataSet, axes3d, indexDistance);
        {
          BinnedED h3d = DistBuilder::Build(evIt->first, 3, pdfConfig, dataSet);
          dist = DistBuilder::Build(evIt->first, 1, pdfConfig, dataSet);
          osc->SetSource(h3d, dist);
        }
        std::cout << "  " << osc->GetName() << ": N_R = " << osc->GetNR() << ", non-zero 3D bins = " << osc->GetNNonZero()
                  << " / " << osc->GetN3DBins() << " (" << 100.0 * osc->GetNNonZero() / osc->GetN3DBins() << "%), (reactor, nu_energy) pairs = "
                  << osc->GetNPairs() << ", events clamped into the nu_energy edge bins = " << nClamped << std::endl;
      }
      else if (evIt->second.GetFlat())
      {
        dist = DistBuilder::BuildFlatDist(evIt->first, num_dimensions, pdfConfig);
      }
      else
      {
        dist = DistBuilder::Build(evIt->first, num_dimensions, pdfConfig, dataSet);
      }
      delete dataSet;
      dist.AddPadding();

      // Save the generated number of events for Beeston Barlow
      genRates[dsIt->first].push_back(dist.Integral());

      // The PDF carries the shape only; its size comes from the rate parameter
      if (dist.Integral() > 0)
        dist.Normalise();
      pdfMap[dsIt->first].push_back(dist);
      BinnedED unscaledPDF = dist;

      normFittingStatuses[dsIt->first]->push_back(INDIRECT);

      // Apply systematics at nominal (Asimov) and fake-data values, group by group as the fit does
      std::vector<double> groupIntegrals;
      dist = ApplyGroups(dist, evIt->second.GetGroup(), systMap[dsIt->first], systGroup, noms, &groupIntegrals);
      BinnedED fakeDataDist = unscaledPDF;
      if (isFakeData)
        fakeDataDist = ApplyGroups(fakeDataDist, evIt->second.GetGroup(), systMap[dsIt->first], systGroup, fdValues);

      std::cout << "  PDF integral after each group:";
      for (size_t iGroup = 0; iGroup < groupIntegrals.size(); iGroup++)
        std::cout << " " << evIt->second.GetGroup().at(iGroup) << " = " << groupIntegrals.at(iGroup);
      std::cout << std::endl;

      // Now scale the Asimov component by expected count, and also save pdf as a histo
      dist.Scale(noms[evIt->first]);
      fakeDataDist.Scale(fdValues[evIt->first]);
      if (dist.GetNDims() != asimov.GetNDims())
      {
        BinnedED marginalised = dist.Marginalise(dataObs);
        asimov.Add(marginalised);
        BinnedED marginalisedPDF = unscaledPDF.Marginalise(dataObs);
        IO::SaveHistogram(marginalised.GetHistogram(), asimovDistDir + "/" + evIt->first + "_" + dsIt->first + ".root", dist.GetName());
        IO::SaveHistogram(marginalisedPDF.GetHistogram(), pdfDir + "/" + evIt->first + "_" + dsIt->first + ".root", unscaledPDF.GetName());
        if (isFakeData)
        {
          BinnedED marginalisedFakeData = fakeDataDist.Marginalise(dataObs);
          fakeDataset.Add(marginalisedFakeData);
          IO::SaveHistogram(marginalisedFakeData.GetHistogram(), fakedataDistDir + "/" + evIt->first + "_" + dsIt->first + ".root", fakeDataDist.GetName());
        }
      }
      else
      {
        asimov.Add(dist);
        IO::SaveHistogram(dist.GetHistogram(), asimovDistDir + "/" + evIt->first + "_" + dsIt->first + ".root", dist.GetName());
        IO::SaveHistogram(unscaledPDF.GetHistogram(), pdfDir + "/" + evIt->first + "_" + dsIt->first + ".root", unscaledPDF.GetName());
        if (isFakeData)
        {
          fakeDataset.Add(fakeDataDist);
          IO::SaveHistogram(fakeDataDist.GetHistogram(), fakedataDistDir + "/" + evIt->first + "_" + dsIt->first + ".root", fakeDataDist.GetName());
        }
      }
    } // End loop over PDFs
    std::cout << std::endl;

    // And save combined histogram (final asimov dataset). This is with everything (including oscillation parameters) nominal.
    IO::SaveHistogram(asimov.GetHistogram(), outDir + "/asimov_" + dsIt->first + ".root", "asimov");
    if (isFakeData)
      IO::SaveHistogram(fakeDataset.GetHistogram(), outDir + "/fakedata_" + dsIt->first + ".root", "fakedata");

    // Now let's load up the data
    BinnedED dataDist;
    if (isAsimov)
    {
      dataDist = asimov;
    }
    else if (isFakeData)
    {
      dataDist = fakeDataset;
    }
    else
    {
      TFile *dataFile = TFile::Open(dataPath[dsIt->first].c_str(), "READ");
      if (!dataFile || dataFile->IsZombie())
        throw std::runtime_error("ERROR: can't open data file " + dataPath[dsIt->first]);

      // Loop through all objects in the file
      TIter nextKey(dataFile->GetListOfKeys());
      TKey *key;
      while ((key = (TKey *)nextKey()))
      {
        std::string className = key->GetClassName();
        std::string objectName = key->GetName();

        if (className == "TNtuple")
        {
          // Load up the data set, and bin the data inside
          ROOTNtuple dataToFit(dataPath[dsIt->first], objectName.c_str());
          dataDist = DistBuilder::Build("data", pdfConfig.GetDataAxisCount(), pdfConfig, (DataSet *)&dataToFit);
          break;
        }
        else if (className == "TH1D")
        {
          TH1D *dataHist = (TH1D *)dataFile->Get(objectName.c_str());
          Histogram loaded = DistTools::ToHist(*dataHist);
          dataDist = BinnedED("data", loaded);
          dataDist.SetObservables(pdfConfig.GetDataBranchNames());
          dataDist.SetAxes(DistBuilder::BuildAxes(pdfConfig, pdfConfig.GetDataAxisCount()));
          break;
        }
      }
      dataFile->Close();
    }
    dataDists[dsIt->first] = dataDist;

    // Now build the likelihood
    BinnedNLLH &lh = testStats.emplace_back();
    lh.SetBuffer("energy", 8, 20);
    lh.SetDataDist(dataDist);
    lh.SetBarlowBeeston(beestonBarlowFlag);
    for (std::map<std::string, Systematic *>::iterator systIt = systMap[dsIt->first].begin(); systIt != systMap[dsIt->first].end(); ++systIt)
      lh.AddSystematic(systIt->second, systGroup[systIt->first]);
    lh.AddPdfs(pdfMap[dsIt->first], pdfGroups[dsIt->first], genRates[dsIt->first], normFittingStatuses[dsIt->first]);

    // Initialise to nominal values
    for (ParameterDict::iterator parIt = mins.begin(); parIt != mins.end(); ++parIt)
    {
      if (std::find(datasetPars[parIt->first].begin(), datasetPars[parIt->first].end(), dsIt->first) == datasetPars[parIt->first].end())
        continue;
      parameterValues[dsIt->first][parIt->first] = noms[parIt->first];
    }

    lh.RegisterFitComponents();
    lh.SetParameters(parameterValues[dsIt->first]);

    std::cout << "Made LLH for Dataset: " << dsIt->first << std::endl
              << std::endl;
  } // End loop over datasets

  // Now combine the LLH for each dataset
  ParameterDict allParVals;
  for (std::map<std::string, ParameterDict>::iterator parMapIt = parameterValues.begin(); parMapIt != parameterValues.end(); parMapIt++)
    for (ParameterDict::iterator parIt = parMapIt->second.begin(); parIt != parMapIt->second.end(); parIt++)
      allParVals.insert(*parIt);

  std::vector<TestStatistic *> rawllhptrs;
  for (BinnedNLLH &lh : testStats)
    rawllhptrs.push_back(&lh);
  StatisticSum fullLLH = Sum(rawllhptrs);

  // Add constraints
  std::vector<std::string> corrPairs;
  for (ParameterDict::iterator corrIt = constrCorrs.begin(); corrIt != constrCorrs.end(); ++corrIt)
  {
    fullLLH.SetConstraint(corrIt->first, constrMeans.at(corrIt->first), constrSigmas.at(corrIt->first), constrCorrParName.at(corrIt->first),
                          constrMeans.at(constrCorrParName.at(corrIt->first)), constrSigmas.at(constrCorrParName.at(corrIt->first)), corrIt->second);
    corrPairs.push_back(constrCorrParName.at(corrIt->first));
  }
  for (ParameterDict::iterator constrIt = constrMeans.begin(); constrIt != constrMeans.end(); ++constrIt)
  {
    // Only add single parameter constraint if correlation hasn't already been applied
    if (constrCorrs.find(constrIt->first) == constrCorrs.end() && std::find(corrPairs.begin(), corrPairs.end(), constrIt->first) == corrPairs.end())
      fullLLH.SetConstraint(constrIt->first, constrIt->second, constrSigmas.at(constrIt->first));
  }
  for (ParameterDict::iterator ratioIt = constrRatioMeans.begin(); ratioIt != constrRatioMeans.end(); ++ratioIt)
    fullLLH.SetConstraint(ratioIt->first, constrRatioParName.at(ratioIt->first), ratioIt->second, constrRatioSigmas.at(ratioIt->first));

  fullLLH.RegisterFitComponents();

  // The MCMC needs a range for every fit parameter and nothing else
  std::set<std::string> fitParNames = fullLLH.GetParameterNames();
  for (const std::string &name : fitParNames)
    if (!mins.count(name) || !maxs.count(name) || !nbins.count(name))
      throw std::runtime_error("ERROR: fit parameter " + name + " needs min, max and nbins in the fit config");

  // Starting point:
  //  - a fixed parameter stays nominal (the sampler never moves it);
  //  - the oscillation parameters (those of the FusedOsc systematics), and any free parameter with no Gaussian
  //    prior on itself, start uniformly in [min, max];
  //  - every other free parameter starts at nominal + Gaussian(0, prior sigma), redrawn until inside [min, max].
  std::set<std::string> oscParNames;
  for (std::map<std::string, std::string>::iterator systIt = systType.begin(); systIt != systType.end(); ++systIt)
    if (IsFusedOsc(systIt->second))
      oscParNames.insert(systParamNames[systIt->first].begin(), systParamNames[systIt->first].end());
  ParameterDict startVals = allParVals;
  std::cout << "Starting point:" << std::endl;
  for (ParameterDict::iterator parIt = startVals.begin(); parIt != startVals.end(); ++parIt)
  {
    const std::string &name = parIt->first;
    std::string rule;
    if (fixedPars.count(name) && fixedPars.at(name))
      rule = "fixed: nominal";
    else if (oscParNames.count(name) || !constrSigmas.count(name))
    {
      parIt->second = Rand::UniformRange(mins.at(name), maxs.at(name));
      rule = oscParNames.count(name) ? "oscillation: uniform in [min, max]" : "free, no Gaussian prior: uniform in [min, max]";
    }
    else
    {
      for (int iDraw = 0; iDraw < 10000; iDraw++)
      {
        const double draw = noms.at(name) + Rand::Gaus(0, constrSigmas.at(name));
        if (draw >= mins.at(name) && draw <= maxs.at(name))
        {
          parIt->second = draw;
          break;
        }
      }
      rule = "nominal + Gaussian(0, prior sigma)";
    }
    std::cout << "  " << name << " = " << parIt->second << "  [" << rule << "]" << std::endl;
  }

  fullLLH.SetParameters(startVals);
  TStopwatch timer;
  timer.Start(true);
  double startLLH = fullLLH.Evaluate();
  timer.Stop();
  std::cout << "NLLH at the starting point: " << startLLH << ", eval time: " << timer.RealTime() << " s" << std::endl;

  // And now the optimiser
  MetropolisSampler sampler;
  for (ParameterDict::iterator it = sigmas.begin(); it != sigmas.end(); ++it)
    sigmas[it->first] *= sigmaScale;
  sampler.SetSigmas(sigmas);
  for (std::map<std::string, bool>::iterator fixIt = fixedPars.begin(); fixIt != fixedPars.end(); ++fixIt)
    if (fixIt->second)
      sampler.Fix(fixIt->first);

  MCMC mh(sampler);
  mh.SetSaveChain(true);
  mh.SetMaxIter(nsteps);
  mh.SetBurnIn(burnin);
  mh.SetMinima(mins);
  mh.SetMaxima(maxs);
  mh.SetInitialTrial(startVals);

  // We're going to minmise not maximise -log(lh) and the
  // mc chain needs to know that the test stat is logged
  // otherwise it will give us the distribution of the log(lh) not the lh
  mh.SetTestStatLogged(true);
  mh.SetFlipSign(true);

  // Create some axes for the mc to fill
  AxisCollection lhAxes;
  for (ParameterDict::iterator it = mins.begin(); it != mins.end(); ++it)
    lhAxes.AddAxis(BinAxis(it->first, it->second, maxs[it->first], nbins[it->first]));
  mh.SetHistogramAxes(lhAxes);

  TStopwatch mcmcTimer;
  mcmcTimer.Start(true);
  const FitResult &res = mh.Optimise(&fullLLH);
  mcmcTimer.Stop();
  std::cout << "MCMC time: " << mcmcTimer.RealTime() << " s for " << nsteps << " steps (" << mcmcTimer.RealTime() / std::max(nsteps, 1) << " s/step)" << std::endl;
  for (std::map<std::string, std::map<std::string, OscWeight1D *>>::iterator dsIt = fusedOscMap.begin(); dsIt != fusedOscMap.end(); ++dsIt)
    for (std::map<std::string, OscWeight1D *>::iterator oscIt = dsIt->second.begin(); oscIt != dsIt->second.end(); ++oscIt)
      std::cout << "  " << oscIt->first << " (" << dsIt->first << ") Construct: " << 1e3 * oscIt->second->GetConstructSeconds() / std::max<size_t>(oscIt->second->GetNConstruct(), 1)
                << " ms/call over " << oscIt->second->GetNConstruct() << " calls" << std::endl;

  ParameterDict bestFit = res.GetBestFit();
  fullLLH.SetParameters(bestFit);

  // Now save the results
  res.SaveAs(outDir + "/fit_result.txt");
  std::cout << "Saved fit result to " << outDir + "/fit_result.txt" << std::endl;

  MCMCSamples samples = mh.GetSamples();

  // Output filename is aboveDirectory_resultDirectory.root, inside outDir
  std::string chainFileName = outDir;
  std::string tempString1 = outDir;
  std::string tempString2 = outDir;
  size_t last_slash_idx = tempString1.find_last_of("/");
  tempString2.erase(0, last_slash_idx + 1);
  if (std::string::npos != last_slash_idx)
  {
    tempString1.erase(last_slash_idx, std::string::npos);
    last_slash_idx = tempString1.find_last_of("/");
    if (std::string::npos != last_slash_idx)
      tempString1.erase(0, last_slash_idx);
  }
  chainFileName = outDir + tempString1 + "_" + tempString2 + ".root";
  TFile *f = new TFile(chainFileName.c_str(), "recreate");
  TTree *outchain = samples.GetChain();
  outchain->Write();
  delete f;

  // Save the histograms, and print the posterior mean and standard deviation of each free parameter
  typedef std::map<std::string, Histogram> HistMap;
  const HistMap &proj1D = samples.Get1DProjections();
  const HistMap &proj2D = samples.Get2DProjections();

  std::cout << "Saving LH projections to \n\t" << projDir1D << "\n\t" << projDir2D << std::endl;
  std::cout << "Posterior (after burn-in, from the 1D projections): mean +- std" << std::endl;
  for (HistMap::const_iterator it = proj1D.begin(); it != proj1D.end(); ++it)
  {
    IO::SaveHistogram(it->second, projDir1D + "/" + it->first + ".root");
    if (fixedPars[it->first])
      continue;
    const BinAxis &axis = it->second.GetAxes().GetAxis(0);
    double sum = 0, sumX = 0, sumXX = 0;
    for (size_t iBin = 0; iBin < axis.GetNBins(); iBin++)
    {
      const double w = it->second.GetBinContent(iBin);
      sum += w;
      sumX += w * axis.GetBinCentre(iBin);
      sumXX += w * axis.GetBinCentre(iBin) * axis.GetBinCentre(iBin);
    }
    if (sum > 0)
      std::cout << "  " << it->first << ": " << sumX / sum << " +- " << sqrt(std::max(0.0, sumXX / sum - sumX * sumX / sum / sum))
                << "   (nominal " << noms[it->first] << ")" << std::endl;
  }
  for (HistMap::const_iterator it = proj2D.begin(); it != proj2D.end(); ++it)
    IO::SaveHistogram(it->second, projDir2D + "/" + it->first + ".root");

  // Now save the postfit distributions for each dataset, applying the best-fit systematics in group order
  std::cout << "Saving scaled histograms and data to \n\t" << postfitDistDir << std::endl;
  for (DSMap::iterator dsIt = dsPDFMap.begin(); dsIt != dsPDFMap.end(); ++dsIt)
  {
    BinnedED postfitDist = BinnedED("postfit dist", systAxes);
    postfitDist.SetObservables(dataObs);
    for (size_t i = 0; i < pdfMap[dsIt->first].size(); i++)
    {
      std::string name = pdfMap[dsIt->first].at(i).GetName();
      BinnedED postfitPDF = ApplyGroups(pdfMap[dsIt->first].at(i), pdfGroups[dsIt->first].at(i), systMap[dsIt->first], systGroup, bestFit);
      postfitPDF.Scale(bestFit[name]);
      if (postfitPDF.GetNDims() != postfitDist.GetNDims())
        postfitPDF = postfitPDF.Marginalise(dataObs);
      IO::SaveHistogram(postfitPDF.GetHistogram(), postfitDistDir + "/" + name + "_" + dsIt->first + ".root");
      // Sum all scaled distributions to get full postfit "dataset"
      postfitDist.Add(postfitPDF);
    }
    IO::SaveHistogram(postfitDist.GetHistogram(), postfitDistDir + "/postfitdist_" + dsIt->first + ".root");

    // And also save the data
    BinnedED dataDist = dataDists[dsIt->first];
    if (dataDist.GetNDims() != postfitDist.GetNDims())
      dataDist = dataDist.Marginalise(dataObs);
    IO::SaveHistogram(dataDist.GetHistogram(), outDir + "/data_" + dsIt->first + ".root");
  }

  // Save autocorrelations
  std::ofstream cofs((outDir + "/auto_correlations.txt").c_str());
  std::vector<double> autocors = samples.GetAutoCorrelations();
  for (size_t i = 0; i < autocors.size(); i++)
    cofs << i << "\t" << autocors.at(i) << "\n";
  cofs.close();

  std::cout << "Fit complete" << std::endl;
}

int main(int argc, char *argv[])
{
  if (argc != 6)
  {
    std::cout << "\nUsage: mcmc <fit_config_file> <eve_config_file> <pdf_config_file> <syst_config_file> <oscgrid_config_file>" << std::endl;
    return 1;
  }

  std::string fitConfigFile(argv[1]);
  std::string eveConfigFile(argv[2]);
  std::string pdfConfigPath(argv[3]);
  std::string systConfigFile(argv[4]);
  std::string oscgridConfigFile(argv[5]);

  mcmc(fitConfigFile, eveConfigFile, pdfConfigPath, systConfigFile, oscgridConfigFile);

  return 0;
}
