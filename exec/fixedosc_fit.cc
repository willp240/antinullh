// Antinu headers
#include <PDFConfigLoader.hh>
#include <DistBuilder.hh>
#include <DistTools.h>
#include <FitConfigLoader.hh>
#include <EventConfigLoader.hh>
#include <SystConfigLoader.hh>
#include <SystFactory.hh>
#include <OscGridConfigLoader.hh>
#include <Utilities.hh>

// OXO headers
#include <ROOTNtuple.h>
#include <BinnedNLLH.h>
#include <StatisticSum.h>
#include <IO.h>
#include <Rand.h>
#include <Minuit.h>

// ROOT headers
#include <TH1D.h>
#include <TKey.h>
#include <TMatrixD.h>
#include <TStopwatch.h>

// c++ headers
#include <sys/stat.h>

using namespace antinufit;

void fixedosc_fit(const std::string &fitConfigFile_,
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
  bool saveOutputs = fitConfig.GetSaveOutputs();
  double tolerance = fitConfig.GetMinuitTolerance();
  int strategy = fitConfig.GetMinuitStrategy();
  std::string method = fitConfig.GetMinuitMethod();
  int iterations = fitConfig.GetIterations();
  std::string outDir = fitConfig.GetOutDir();
  ParameterDict constrMeans = fitConfig.GetConstrMeans();
  ParameterDict constrSigmas = fitConfig.GetConstrSigmas();
  ParameterDict mins = fitConfig.GetMinima();
  ParameterDict maxs = fitConfig.GetMaxima();
  ParameterDict noms = fitConfig.GetNominals();
  ParameterDict sigmas = fitConfig.GetSigmas();
  ParameterDict constrRatioMeans = fitConfig.GetConstrRatioMeans();
  ParameterDict constrRatioSigmas = fitConfig.GetConstrRatioSigmas();
  std::map<std::string, std::string> constrRatioParName = fitConfig.GetConstrRatioParName();
  ParameterDict constrCorrs = fitConfig.GetConstrCorrs();
  std::map<std::string, std::string> constrCorrParName = fitConfig.GetConstrCorrParName();
  ParameterDict fdValues = fitConfig.GetFakeDataVals();

  std::string pdfDir = outDir + "/unscaled_pdfs";
  std::string asimovDistDir = outDir + "/asimov_dists";
  std::string fakedataDistDir = outDir + "/fakedata_dists";
  std::string postfitDistDir = outDir + "/postfit_dists";
  if (saveOutputs)
  {
    // Create output directories
    struct stat st = {0};
    if (stat(outDir.c_str(), &st) == -1)
      mkdir(outDir.c_str(), 0700);
    if (stat(pdfDir.c_str(), &st) == -1)
      mkdir(pdfDir.c_str(), 0700);
    if (stat(asimovDistDir.c_str(), &st) == -1)
      mkdir(asimovDistDir.c_str(), 0700);
    if (stat(fakedataDistDir.c_str(), &st) == -1)
      mkdir(fakedataDistDir.c_str(), 0700);
    if (stat(postfitDistDir.c_str(), &st) == -1)
      mkdir(postfitDistDir.c_str(), 0700);
  }

  // Load up all the event types we want to contribute
  typedef std::map<std::string, EventConfig> EvMap;
  typedef std::map<std::string, std::map<std::string, EventConfig>> DSMap;
  EventConfigLoader evLoader(evConfigFile_);
  DSMap dsPDFMap = evLoader.LoadActive();
  std::map<std::string, std::string> dataPath = evLoader.GetDataPaths();
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

  // Load up the oscillation probability grids
  OscGridConfigLoader oscGridLoader(oscGridConfigFile_);
  OscGridConfig oscGridConfig = oscGridLoader.Load();
  std::map<int, OscGrid *> dummyOscGridMap;

  // First read the reactor distance info
  std::string reactorjson = oscGridConfig.GetReactorsJsonFile();
  std::unordered_map<int, double> indexDistance = LoadIndexDistanceMap(reactorjson);

  // Loop over systematics and declare each one
  std::map<std::string, std::map<std::string, Systematic *>> systMap;
  for (std::map<std::string, std::string>::iterator systIt = systType.begin(); systIt != systType.end(); ++systIt)
  {
    // Check any parameter defined for a systematic, was declared in the fit config
    for (std::map<std::string, std::vector<std::string>>::iterator paramIt = systParamNames.begin(); paramIt != systParamNames.end(); ++paramIt)
    {
      for (int iParam = 0; iParam < paramIt->second.size(); iParam++)
      {
        if (noms.find(paramIt->second.at(iParam)) == noms.end())
        {
          std::cout << "ERROR: Syst config defines systematic with parameter " << paramIt->second.at(iParam) << ", but this parameter is not defined in the fit config" << std::endl;
          throw;
        }
      }
    }
    fullParamNameVec.insert(fullParamNameVec.end(), systParamNames[systIt->first].begin(), systParamNames[systIt->first].end());

    // Now build the systematic
    Systematic *syst = SystFactory::New(systIt->first, systType[systIt->first], systParamNames[systIt->first], noms, dummyOscGridMap, indexDistance);
    AxisCollection systAxes = DistBuilder::BuildAxes(pdfConfig, systDistObs[systIt->first].size());
    syst->SetAxes(systAxes);
    // The "dimensions" the systematic applies to
    syst->SetTransformationObs(systTransObs[systIt->first]);
    // All the "dimensions" of the dataset
    syst->SetDistributionObs(systDistObs[systIt->first]);
    syst->Construct();
    for (std::vector<std::string>::iterator dsIt = systDataSets[systIt->first].begin(); dsIt != systDataSets[systIt->first].end(); ++dsIt)
      systMap[*dsIt][systIt->first] = syst;
  }

  // Check we've got the oscillation parameters we need, and whether theta is sin-ed or not
  bool hasTheta12 = noms.find("theta12") != noms.end();
  bool hasSinTheta12 = noms.find("sintheta12") != noms.end();
  bool hasSinSqTheta12 = noms.find("sinsqtheta12") != noms.end();
  bool hasDeltam21 = noms.find("deltam21") != noms.end();
  int thetacount = int(hasTheta12) + int(hasSinTheta12) + int(hasSinSqTheta12);
  if (thetacount == 0 || !hasDeltam21)
  {
    throw std::runtime_error("ERROR: A theta12, sintheta12, or sinsqtheta12 parameter, along with a deltam21 parameter, must be provided in the fitconfig.");
  }
  if (thetacount > 1)
  {
    throw std::runtime_error("ERROR: More than one of theta12, sintheta12, sinsqtheta12 parameters were set in the fitconfig. There must be exactly one.");
  }
  std::string theta12name;
  if (hasTheta12)
    theta12name = "theta12";
  else if (hasSinTheta12)
    theta12name = "sintheta12";
  else if (hasSinSqTheta12)
    theta12name = "sinsqtheta12";

  double deltam21 = noms["deltam21"];
  double theta12 = noms[theta12name];
  std::cout << std::endl;

  // A parameter could have been defined in the fit config but isn't associated with a pdf or systematic
  // If so ignore that parameter
  ParameterDict::iterator parIt = noms.begin();
  size_t numParams = noms.size();

  std::map<std::string, std::vector<std::string>> datasetPars;
  for (int iParam = 0; iParam < numParams; iParam++)
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

  PrintParams(mins, maxs, noms, constrMeans, constrSigmas, constrRatioMeans, constrRatioSigmas, constrRatioParName, constrCorrs, constrCorrParName, datasetPars);

  // Create the individual PDFs and Asimov components, for each dataset, and make the component LLH objects
  std::map<std::string, std::vector<BinnedED>> pdfMap;
  std::map<std::string, std::vector<std::vector<std::string>>> pdfGroups;
  std::vector<BinnedNLLH> testStats;
  std::map<std::string, ParameterDict> parameterValues;
  std::map<std::string, double > reactorRatio; // Ratio of oscillated to unoscillated number of reactor IBDs
  std::map<std::string, double > reactorRatioFD; // Ratio of oscillated to unoscillated number of reactor IBDs for the fake dataset if we're making one
  std::map<std::string, std::vector<int>> genRates;
  std::map<std::string, std::vector<NormFittingStatus> *> normFittingStatuses;

  // Now we're going to loop over datasets and build the asimov datsets
  for (DSMap::iterator dsIt = dsPDFMap.begin(); dsIt != dsPDFMap.end(); ++dsIt)
  {
    std::cout << "Building Asimov for dataset: " << dsIt->first << std::endl;
    std::vector<NormFittingStatus> *normStatVec = new std::vector<NormFittingStatus>;
    normFittingStatuses[dsIt->first] = normStatVec;

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
      pdfGroups[dsIt->first].push_back(evIt->second.GetGroup());

      // Get MC events from file
      DataSet *dataSet;
      try
      {
        std::cout << "Loading " << evIt->second.GetPrunedPath() << std::endl;
        dataSet = new ROOTNtuple(evIt->second.GetPrunedPath(), "pruned");
      }
      catch (const IOError &e_)
      {
        std::cout << "Warning: skipping " << evIt->first << " couldn't open files:\n\t" << e_.what() << std::endl;
        continue;
      }

      // Build distribution of those events
      BinnedED dist;
      BinnedED fakeDataDist;
      BinnedED unscaledPDF;
      int num_dimensions = evIt->second.GetNumDimensions();

      if (evIt->first.find("reactor_nubar") != std::string::npos)
      {
        reactorRatio[evIt->first];
        reactorRatioFD[evIt->first];
        // Whichever form of theta12 we have, let's make a sin^2(theta12) to hand to the oscillated dist builder
        double theta12_param = theta12;
        double theta12_fdparam = fdValues[theta12name];
        if (hasSinTheta12)
        {
          theta12_param = theta12 * theta12;
          theta12_fdparam = fdValues[theta12name] * fdValues[theta12name];
        }
        else if (hasTheta12)
        {
          theta12_param = sin(M_PI * theta12 / 180) * sin(M_PI * theta12 / 180);
          theta12_fdparam = sin(M_PI * fdValues[theta12name] / 180) * sin(M_PI * fdValues[theta12name] / 180);
        }

        // Build the distribution with oscillation parameters at their nominal values
        dist = DistBuilder::BuildOscillatedDist(evIt->first, num_dimensions, pdfConfig, dataSet, deltam21, theta12_param, indexDistance, reactorRatio[evIt->first]);
        fakeDataDist = DistBuilder::BuildOscillatedDist(evIt->first, num_dimensions, pdfConfig, dataSet, fdValues["deltam21"], theta12_fdparam, indexDistance, reactorRatioFD[evIt->first]);

        // Now we will scale the constraint on the unoscillated reactor flux by the ratio of the oscillated to unoscillated number of events
        if (constrMeans.find(evIt->first) != constrMeans.end())
        {
          constrMeans[evIt->first] = constrMeans[evIt->first] * reactorRatio[evIt->first];
          constrSigmas[evIt->first] = constrSigmas[evIt->first] * reactorRatio[evIt->first];
        }
        noms[evIt->first] = noms[evIt->first] * reactorRatio[evIt->first];
        mins[evIt->first] = mins[evIt->first] * reactorRatio[evIt->first];
        maxs[evIt->first] = maxs[evIt->first] * reactorRatio[evIt->first];

        fdValues[evIt->first] = fdValues[evIt->first] * reactorRatioFD[evIt->first];
      }
      else
      {
        // For all other PDFs, just use Build
        dist = DistBuilder::Build(evIt->first, num_dimensions, pdfConfig, dataSet);
        fakeDataDist = DistBuilder::Build(evIt->first, num_dimensions, pdfConfig, dataSet);
      }
      dist.AddPadding();
      fakeDataDist.AddPadding();

      // Save the generated number of events for Beeston Barlow
      genRates[dsIt->first].push_back(dist.Integral());

      // Scale for PDF and add to vector
      if (dist.Integral() && fakeDataDist.Integral())
      {
        dist.Normalise();
        fakeDataDist.Normalise();
      }

      pdfMap[dsIt->first].push_back(dist);
      unscaledPDF = dist;

      normFittingStatuses[dsIt->first]->push_back(INDIRECT);

      // Apply nominal systematic variables
      for (std::map<std::string, Systematic *>::iterator systIt = systMap[dsIt->first].begin(); systIt != systMap[dsIt->first].end(); ++systIt)
      {
        // If group is "", we apply to all groups
        if (systGroup[systIt->first] == "" || std::find(pdfGroups[dsIt->first].back().begin(), pdfGroups[dsIt->first].back().end(), systGroup[systIt->first]) != pdfGroups[dsIt->first].back().end())
        {
          double norm;
          dist = systIt->second->operator()(dist, &norm);

          // Set syst parameter to fake data value, and apply to fake data dist and rescale
          std::set<std::string> systParamNames = systIt->second->GetParameterNames();
          for (auto itSystParam = systParamNames.begin(); itSystParam != systParamNames.end(); ++itSystParam)
          {
            systIt->second->SetParameter(*itSystParam, fdValues[*itSystParam]);
          }
          systIt->second->Construct();
          fakeDataDist = systIt->second->operator()(fakeDataDist, &norm);
        }
      }

      // Now scale the Asimov component by expected count, and also save pdf as a histo
      dist.Scale(noms[evIt->first]);
      fakeDataDist.Scale(fdValues[evIt->first]);
      if (dist.GetNDims() != asimov.GetNDims())
      {
        BinnedED marginalised = dist.Marginalise(dataObs);
        asimov.Add(marginalised);
        BinnedED marginalisedPDF = unscaledPDF.Marginalise(dataObs);
        if (saveOutputs)
        {
          IO::SaveHistogram(marginalised.GetHistogram(), asimovDistDir + "/" + evIt->first + "_" + dsIt->first + ".root", dist.GetName());
          IO::SaveHistogram(marginalisedPDF.GetHistogram(), pdfDir + "/" + evIt->first + "_" + dsIt->first + ".root", unscaledPDF.GetName());
        }
        // Also scale fake data dist by fake data value
        if (isFakeData)
        {
          BinnedED marginalisedFakeData = fakeDataDist.Marginalise(dataObs);
          fakeDataset.Add(marginalisedFakeData);
          if (saveOutputs)
            IO::SaveHistogram(marginalisedFakeData.GetHistogram(), fakedataDistDir + "/" + evIt->first + "_" + dsIt->first + ".root", fakeDataDist.GetName());
        }
      }
      else
      {
        asimov.Add(dist);
        if (saveOutputs)
        {
          IO::SaveHistogram(dist.GetHistogram(), asimovDistDir + "/" + evIt->first + "_" + dsIt->first + ".root", dist.GetName());
          IO::SaveHistogram(unscaledPDF.GetHistogram(), pdfDir + "/" + evIt->first + "_" + dsIt->first + ".root", unscaledPDF.GetName());
        }
        // Also scale fake data dist by fake data value
        if (isFakeData)
        {
          fakeDataset.Add(fakeDataDist);
          if (saveOutputs)
            IO::SaveHistogram(fakeDataDist.GetHistogram(), fakedataDistDir + "/" + evIt->first + "_" + dsIt->first + ".root", fakeDataDist.GetName());
        }
      }
    } // End loop over PDFs
    std::cout << std::endl;

    // And save combined histogram (final asimov dataset). This is with everything (including oscillation parameters) nominal. This is what we
    // compare to to calculate the llh
    if (saveOutputs)
    {
      IO::SaveHistogram(asimov.GetHistogram(), outDir + "/asimov_" + dsIt->first + ".root", "asimov");
      if (isFakeData)
      {
        IO::SaveHistogram(fakeDataset.GetHistogram(), outDir + "/fakedata_" + dsIt->first + ".root", "fakedata");
      }
    }
    std::cout << std::endl;

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

      // Could be h5 or root file
      if (dataPath[dsIt->first].substr(dataPath[dsIt->first].find_last_of(".") + 1) == "h5")
      {
        Histogram loaded = IO::LoadHistogram(dataPath[dsIt->first]);
        dataDist = BinnedED("data", loaded);
        dataDist.SetObservables(pdfConfig.GetBranchNames());
      }
      else // If a root file
      {

        TFile *dataFile = TFile::Open(dataPath[dsIt->first].c_str(), "READ");
        if (!dataFile || dataFile->IsZombie())
        {
          std::cerr << "Error opening data file " << dataFile << std::endl;
          throw;
        }

        // Get the list of keys in the file
        TList *keyList = dataFile->GetListOfKeys();
        if (!keyList)
        {
          std::cerr << "Error: No keys found in the datafile " << dataPath[dsIt->first] << std::endl;
          throw;
        }

        // Loop through all objects in the file
        TIter nextKey(keyList);
        TKey *key;
        while ((key = (TKey *)nextKey()))
        {
          std::string className = key->GetClassName();
          std::string objectName = key->GetName();

          if (className == "TNtuple")
          {
            // Load up the data set
            ROOTNtuple dataToFit(dataPath[dsIt->first], objectName.c_str());

            // And bin the data inside
            dataDist = DistBuilder::Build("data", pdfConfig.GetDataAxisCount(), pdfConfig, (DataSet *)&dataToFit);
            break; // Stop once we find an nTuple
          }
          else if (className == "TH1D")
          {
            TH1D *dataHist = (TH1D *)dataFile->Get(objectName.c_str());
            Histogram loaded = DistTools::ToHist(*dataHist);
            dataDist = BinnedED("data", loaded);
            dataDist.SetObservables(pdfConfig.GetDataBranchNames());
            AxisCollection axes = DistBuilder::BuildAxes(pdfConfig, pdfConfig.GetDataAxisCount());
            dataDist.SetAxes(axes);
            break;
          }
        }
      }
    }
    dataDists[dsIt->first] = dataDist;

    // Now build the likelihood
    BinnedNLLH lh;
    lh.SetBuffer("energy", 1, 20);
    // Add our data
    lh.SetDataDist(dataDist);
    // Set whether or not to use Beeston Barlow
    lh.SetBarlowBeeston(beestonBarlowFlag);
    // Add the systematics and any prior constraints
    for (std::map<std::string, Systematic *>::iterator systIt = systMap[dsIt->first].begin(); systIt != systMap[dsIt->first].end(); ++systIt)
      lh.AddSystematic(systIt->second, systGroup[systIt->first]);
    // Add our pdfs
    lh.AddPdfs(pdfMap[dsIt->first], pdfGroups[dsIt->first], genRates[dsIt->first], normFittingStatuses[dsIt->first]);

    // And constraints
    std::vector<std::string> corrPairs;
    for (ParameterDict::iterator corrIt = constrCorrs.begin(); corrIt != constrCorrs.end(); ++corrIt)
    {
      // If either parameter doesn't exist for this dataset, move along
      if (std::find(datasetPars[corrIt->first].begin(), datasetPars[corrIt->first].end(), dsIt->first) == datasetPars[corrIt->first].end() ||
          std::find(datasetPars[constrCorrParName.at(corrIt->first)].begin(), datasetPars[constrCorrParName.at(corrIt->first)].end(), dsIt->first) == datasetPars[constrCorrParName.at(corrIt->first)].end())
        continue;
      lh.SetConstraint(corrIt->first, constrMeans.at(corrIt->first), constrSigmas.at(corrIt->first), constrCorrParName.at(corrIt->first),
                       constrMeans.at(constrCorrParName.at(corrIt->first)), constrSigmas.at(constrCorrParName.at(corrIt->first)), corrIt->second);
      corrPairs.push_back(constrCorrParName.at(corrIt->first));
    }
    for (ParameterDict::iterator constrIt = constrMeans.begin(); constrIt != constrMeans.end(); ++constrIt)
    {
      // If parameter doesn't exist for this dataset, move along
      if (std::find(datasetPars[constrIt->first].begin(), datasetPars[constrIt->first].end(), dsIt->first) == datasetPars[constrIt->first].end())
        continue;

      // Only add single parameter constraint if correlation hasn't already been applied
      if (constrCorrs.find(constrIt->first) == constrCorrs.end() && std::find(corrPairs.begin(), corrPairs.end(), constrIt->first) == corrPairs.end())
        lh.SetConstraint(constrIt->first, constrIt->second, constrSigmas.at(constrIt->first));
    }
    for (ParameterDict::iterator ratioIt = constrRatioMeans.begin(); ratioIt != constrRatioMeans.end(); ++ratioIt)
    {
      // If parameter doesn't exist for this dataset, move along
      if (std::find(datasetPars[ratioIt->first].begin(), datasetPars[ratioIt->first].end(), dsIt->first) == datasetPars[ratioIt->first].end())
        continue;

      lh.SetConstraint(ratioIt->first, constrRatioParName.at(ratioIt->first), ratioIt->second, constrRatioSigmas.at(ratioIt->first));
    }

    // And finally bring it all together
    lh.RegisterFitComponents();

    // Initialise to nominal values
    for (ParameterDict::iterator parIt = mins.begin(); parIt != mins.end(); ++parIt)
    {
      if (std::find(datasetPars[parIt->first].begin(), datasetPars[parIt->first].end(), dsIt->first) == datasetPars[parIt->first].end())
        continue;
      parameterValues[dsIt->first][parIt->first] = noms[parIt->first];
      if (isFakeData)
      {
        parameterValues[dsIt->first][parIt->first] = fdValues[parIt->first];
      }
    }
    // If we have constraints, initialise to those
    for (ParameterDict::iterator constrIt = constrMeans.begin(); constrIt != constrMeans.end(); ++constrIt)
    {
      if (std::find(datasetPars[constrIt->first].begin(), datasetPars[constrIt->first].end(), dsIt->first) == datasetPars[constrIt->first].end())
        continue;
      parameterValues[dsIt->first][constrIt->first] = constrMeans[constrIt->first];
    }
    // Set to these initial values
    lh.SetParameters(parameterValues[dsIt->first]);
    testStats.push_back(std::move(lh));
    std::cout << "Made LLH for Datset: " << dsIt->first << std::endl << std::endl;
  } // End loop over datasets

  // Now combine the LLH for each dataset
  ParameterDict allParVals;
  for (std::map<std::string, ParameterDict>::iterator parMapIt = parameterValues.begin(); parMapIt != parameterValues.end(); parMapIt++)
  {
    for (ParameterDict::iterator parIt = parMapIt->second.begin(); parIt != parMapIt->second.end(); parIt++)
      allParVals[parIt->first] = parIt->second;
  }

  std::vector<TestStatistic *> rawllhptrs;
  for (BinnedNLLH &lh : testStats)
  {
    rawllhptrs.push_back(&lh);
  }
  StatisticSum fullLLH = Sum(rawllhptrs);
  fullLLH.RegisterFitComponents();

  TStopwatch timer;
  timer.Start( true );
  fullLLH.Evaluate();
  timer.Stop();
  std::cout << "Eval time: " << timer.RealTime() << std::endl;

  // Now do a fit!
  Minuit min;
  std::cout << method << " " << iterations << " " << tolerance << " " << strategy << std::endl;
  min.SetMethod(method); // Simplex
  min.SetMaxCalls(iterations); // 1200
  min.SetTolerance(tolerance);
  min.SetStrategy(strategy);
  min.SetMinima(mins);
  min.SetMaxima(maxs);
  min.SetInitialValues(allParVals);
  min.SetInitialErrors(sigmas);
  std::cout << "Run rabbit run!" << std::endl;
  TStopwatch minuitTimer;
  minuitTimer.Start( true );
  FitResult res = min.Optimise(&fullLLH);
  minuitTimer.Stop();
  std::cout << "Minuit time: " << minuitTimer.RealTime() << std::endl;
  res.SetPrintPrecision(4);
  res.Print();
  ParameterDict bestFit = res.GetBestFit();
  fullLLH.SetParameters(bestFit);
  bool validFit = res.GetValid();
  double finalLLH = fullLLH.Evaluate();

  // Now save the results
  if (saveOutputs)
  {
    res.SaveAs(outDir + "/fit_result.txt");
    std::ofstream file(outDir + "/fit_result.txt", std::ios::app);
    file << "\nLLH: " << finalLLH << "\n";
    file << "\nFit Valid: " << validFit << std::endl;
    file.close();
    TFile *outFile = new TFile((outDir + "/fit_result.root").c_str(), "RECREATE");
    DenseMatrix covMatrix = res.GetCovarianceMatrix();
    std::vector<std::string> paramNames;
    std::vector<double> paramVals;
    std::vector<double> paramErr;
    TMatrixD covTMatrixD(bestFit.size(), bestFit.size());
    for (ParameterDict::iterator parIt = bestFit.begin(); parIt != bestFit.end(); ++parIt)
    {
      paramNames.push_back(parIt->first);
      paramVals.push_back(parIt->second);
      if (validFit)
      {
        paramErr.push_back(sqrt(covMatrix.GetComponent(paramNames.size() - 1, paramNames.size() - 1)));
        for (int iParam = 0; iParam < paramNames.size(); iParam++)
        {
          covTMatrixD[paramNames.size() - 1][iParam] = covMatrix.GetComponent(paramNames.size() - 1, iParam);
          covTMatrixD[iParam][paramNames.size() - 1] = covMatrix.GetComponent(iParam, paramNames.size() - 1);
        }
      }
    }
    paramNames.push_back("LLH");
    paramVals.push_back(finalLLH);
    paramNames.push_back("FitValid");
    paramVals.push_back(validFit);
    outFile->WriteObject(&paramNames, "paramNames");
    outFile->WriteObject(&paramVals, "paramVals");
    outFile->WriteObject(&paramErr, "paramErr");
    outFile->WriteObject(&covTMatrixD, "covMatrix");
    std::cout << "Saved fit result to " << outDir + "/fit_result.txt and " << outDir << "/fit_result.root" << std::endl;

    // Now save the postfit distributions for each dataset
    for (DSMap::iterator dsIt = dsPDFMap.begin(); dsIt != dsPDFMap.end(); ++dsIt)
    {
      // Initialise postfit distributions to same axis as data
      BinnedED postfitDist = BinnedED("postfit dist", systAxes);
      postfitDist.SetObservables(dataObs);

      // Scale the distributions to the correct heights. They are named the same as their fit parameters
      std::cout << "Saving scaled histograms and data for " << dsIt->first << " to \n\t" << postfitDistDir << std::endl;

      if (dataDists[dsIt->first].GetHistogram().GetNDims() < 3)
      {
        ParameterDict bestFit = res.GetBestFit();

        for (size_t i = 0; i < pdfMap[dsIt->first].size(); i++)
        {
          std::string name = pdfMap[dsIt->first].at(i).GetName();
          pdfMap[dsIt->first][i].Normalise();
          // Apply bestfit systematic variables
          for (std::map<std::string, Systematic *>::iterator systIt = systMap[dsIt->first].begin(); systIt != systMap[dsIt->first].end(); ++systIt)
          {
            // If group is "", we apply to all groups
            if (systGroup[systIt->first] == "" || std::find(pdfGroups[dsIt->first].back().begin(), pdfGroups[dsIt->first].back().end(), systGroup[systIt->first]) != pdfGroups[dsIt->first].back().end())
            {
              std::set<std::string> systParamNames = systIt->second->GetParameterNames();
              for (auto itSystParam = systParamNames.begin(); itSystParam != systParamNames.end(); ++itSystParam)
              {
                systIt->second->SetParameter(*itSystParam, bestFit[*itSystParam]);
              }
              double norm;
              systIt->second->Construct();
              pdfMap[dsIt->first][i] = systIt->second->operator()(pdfMap[dsIt->first][i], &norm);
            }
          }
          pdfMap[dsIt->first][i].Scale(bestFit[name]);
          IO::SaveHistogram(pdfMap[dsIt->first][i].GetHistogram(), postfitDistDir + "/" + name + "_" + dsIt->first + ".root");
          // Sum all scaled distributions to get full postfit "dataset"
          postfitDist.Add(pdfMap[dsIt->first][i]);
        }
        // WP: name should include dataset name
        IO::SaveHistogram(postfitDist.GetHistogram(), postfitDistDir + "/postfitdist_" + dsIt->first + ".root");
      }
      else
      {
        ParameterDict bestFit = res.GetBestFit();
        for (size_t i = 0; i < pdfMap[dsIt->first].size(); i++)
        {
          std::string name = pdfMap[dsIt->first].at(i).GetName();
          pdfMap[dsIt->first][i].Normalise();
          // Apply bestfit systematic variables
          for (std::map<std::string, Systematic *>::iterator systIt = systMap[dsIt->first].begin(); systIt != systMap[dsIt->first].end(); ++systIt)
          {
            // If group is "", we apply to all groups
            if (systGroup[systIt->first] == "" || std::find(pdfGroups[dsIt->first].back().begin(), pdfGroups[dsIt->first].back().end(), systGroup[systIt->first]) != pdfGroups[dsIt->first].back().end())
            {
              std::set<std::string> systParamNames = systIt->second->GetParameterNames();
              for (auto itSystParam = systParamNames.begin(); itSystParam != systParamNames.end(); ++itSystParam)
              {
                systIt->second->SetParameter(*itSystParam, bestFit[*itSystParam]);
              }
              double norm;
              systIt->second->Construct();
              pdfMap[dsIt->first][i] = systIt->second->operator()(pdfMap[dsIt->first][i], &norm);
            }
          }
          pdfMap[dsIt->first][i].Scale(bestFit[name]);

          std::vector<std::string> keepObs;
          keepObs.push_back("energy");
          pdfMap[dsIt->first][i] = pdfMap[dsIt->first][i].Marginalise(keepObs);
          IO::SaveHistogram(pdfMap[dsIt->first][i].GetHistogram(), postfitDistDir + "/" + name + "_" + dsIt->first + ".root");
          // Sum all scaled distributions to get full postfit "dataset"
          postfitDist.Add(pdfMap[dsIt->first][i]);
        }

        IO::SaveHistogram(postfitDist.GetHistogram(), postfitDistDir + "/postfitdist_" + dsIt->first + ".root");
      }

      // And also save the data
      if (dataDists[dsIt->first].GetHistogram().GetNDims() < 3)
      {
        IO::SaveHistogram(dataDists[dsIt->first].GetHistogram(), outDir + "/" + "data_" + dsIt->first + ".root");
      }
      else
      {
        std::vector<std::string> keepObs;
        keepObs.push_back("energy");
        dataDists[dsIt->first] = dataDists[dsIt->first].Marginalise(keepObs);
        // WP: name will include ds name
        IO::SaveHistogram(dataDists[dsIt->first].GetHistogram(), outDir + "/" + "data_" + dsIt->first + ".root");
      }

    } // End loop over datasets

  } // End if saving outputs

  std::cout << "Fit complete for:" << std::endl;
  std::cout << "deltam: " << deltam21 << std::endl;
  std::cout << theta12name << ": " << theta12 << std::endl;
  std::cout << "LLH: " << finalLLH << std::endl;
  for (auto reacRatioIt = reactorRatio.begin(); reacRatioIt != systParamNames.end(); ++reacRatioIt)
  {
    std::cout << reacRatioIt->first << " Reactor Ratio: " << reacRatioIt->second << std::endl;
  }
  std::cout << "FitValid: " << validFit << std::endl;
  std::cout << std::endl
            << std::endl;
}

int main(int argc, char *argv[])
{
  if (argc != 6)
  {
    std::cout << "\nUsage fixedosc_fit <fit_config_file> <eve_config_file> <pdf_config_file> <syst_config_file> <oscgrid_config_file>" << std::endl;
    return 1;
  }

  std::string fitConfigFile(argv[1]);
  std::string eveConfigFile(argv[2]);
  std::string pdfConfigPath(argv[3]);
  std::string systConfigFile(argv[4]);
  std::string oscgridConfigFile(argv[5]);

  fixedosc_fit(fitConfigFile, eveConfigFile, pdfConfigPath, systConfigFile, oscgridConfigFile);

  return 0;
}
