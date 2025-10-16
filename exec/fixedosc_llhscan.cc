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

// ROOT headers
#include <TH1D.h>
#include <TKey.h>

// c++ headers
#include <sys/stat.h>

using namespace antinufit;

void fixedosc_llhscan(const std::string &fitConfigFile_,
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
  std::string outDir = fitConfig.GetOutDir();
  ParameterDict constrMeans = fitConfig.GetConstrMeans();
  ParameterDict constrSigmas = fitConfig.GetConstrSigmas();
  ParameterDict mins = fitConfig.GetMinima();
  ParameterDict maxs = fitConfig.GetMaxima();
  ParameterDict noms = fitConfig.GetNominals();
  ParameterDict constrRatioMeans = fitConfig.GetConstrRatioMeans();
  ParameterDict constrRatioSigmas = fitConfig.GetConstrRatioSigmas();
  std::map<std::string, std::string> constrRatioParName = fitConfig.GetConstrRatioParName();
  ParameterDict constrCorrs = fitConfig.GetConstrCorrs();
  std::map<std::string, std::string> constrCorrParName = fitConfig.GetConstrCorrParName();
  ParameterDict fdValues = fitConfig.GetFakeDataVals();
  std::map<std::string, bool> fixedPars = fitConfig.GetFixPars();
  std::map<std::string, std::string> labelName = fitConfig.GetTexLabels();

  std::string pdfDir = outDir + "/unscaled_pdfs";
  std::string asimovDistDir = outDir + "/asimov_dists";
  std::string fakedataDistDir = outDir + "/fakedata_dists";
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

  // Load up the oscillation probability grids information
  // This is really just to get the reactors json file
  OscGridConfigLoader oscGridLoader(oscGridConfigFile_);
  OscGridConfig oscGridConfig = oscGridLoader.Load();
  std::map<int, OscGrid *> dummyOscGridMap;

  // Read the reactor distance info
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

  // Store the oscillation parameter scan region. We got these from the fit config file, but are about to delete them as they're
  // not used in the OXO llh
  double deltam21_nom = noms["deltam21"];
  double deltam21_min = mins["deltam21"];
  double deltam21_max = maxs["deltam21"];
  double deltam21_fd = fdValues["deltam21"];
  double theta12_nom = noms[theta12name];
  double theta12_min = mins[theta12name];
  double theta12_max = maxs[theta12name];
  double theta12_fd = fdValues[theta12name];

  std::cout << std::endl;

  // Define the number of points
  int npoints = 150;
  int countwidth = double(npoints) / double(5);

  // A parameter could have been defined in the fit config but isn't associated with a pdf or systematic
  // If so ignore that parameter
  ParameterDict::iterator parIt = noms.begin();
  size_t numParams = noms.size();

  std::map<std::string, int> reacPDFNum;
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
  std::map<std::string, ParameterDict> parameterValues;
  std::map<std::string, std::vector<BinnedED>> oscPDFsMap;
  double reactorRatio = 1.0; // Ratio of oscillated to unoscillated number of reactor IBDs
  std::map<std::string, std::vector<int>> genRates;
  std::map<std::string, std::vector<NormFittingStatus> *> normFittingStatuses;

  // Now we're going to loop over datasets and build the asimov datsets
  for (DSMap::iterator dsIt = dsPDFMap.begin(); dsIt != dsPDFMap.end(); ++dsIt)
  {
    std::cout << "Building Asimov for dataset: " << dsIt->first << std::endl;
    reacPDFNum[dsIt->first] = -999;
    std::vector<NormFittingStatus> *normStatVec = new std::vector<NormFittingStatus>;
    normFittingStatuses[dsIt->first] = normStatVec;

    // Create the empty full dist
    BinnedED asimov = BinnedED("asimov", systAxes);
    asimov.SetObservables(dataObs);
    // And an empty fake data dist
    BinnedED fakeDataset = BinnedED("fake_dataset", systAxes);
    fakeDataset.SetObservables(dataObs);

    // And we're going to create a scaled "asimov" dataset and reactor PDF for each point in the oscillation parameter scans
    // The oscDatasets aren't actually used in the scan, but we save them as they're useful to look at. But the teststat object
    // should build exact equivalents if everything has gone to plan
    std::vector<BinnedED> oscDatasets;
    for (int iOscPoints = 0; iOscPoints < 2 * npoints; iOscPoints++)
    {
      BinnedED oscAsmv = BinnedED("asimov", systAxes);
      oscAsmv.SetObservables(dataObs);
      oscDatasets.push_back(oscAsmv);
    }

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

      // If it's the reactor PDF, two things are different from the others. Firstly we use the BuildOscillatedDist rather than just Build.
      // Also, we're going to loop over oscillation scan points to make the reactor's PDF for that point as well
      if (evIt->first.find("reactor_nubar") != std::string::npos)
      {
        reacPDFNum[dsIt->first] = pdfMap[dsIt->first].size();

        // Whichever form of theta12 we have, let's make a sin^2(theta12) to hand to the oscillated dist builder
        double theta12_param = theta12_nom;
        double theta12_fdparam = theta12_fd;
        if (hasSinTheta12)
        {
          theta12_param = theta12_nom * theta12_nom;
          theta12_fdparam = theta12_fd * theta12_fd;
          mins[theta12name] = mins[theta12name] * mins[theta12name];
          maxs[theta12name] = maxs[theta12name] * maxs[theta12name];
        }
        else if (hasTheta12)
        {
          theta12_param = sin(M_PI * theta12_nom / 180) * sin(M_PI * theta12_nom / 180);
          theta12_fdparam = sin(M_PI * theta12_fd / 180) * sin(M_PI * theta12_fd / 180);
          mins[theta12name] = sin(M_PI * mins[theta12name] / 180) * sin(M_PI * mins[theta12name] / 180);
          maxs[theta12name] = sin(M_PI * maxs[theta12name] / 180) * sin(M_PI * maxs[theta12name] / 180);
        }

        // Build the distribution with oscillation parameters at their nominal values
        // Pass double by reference here to get ratio of unoscillated to oscillated events
        dist = DistBuilder::BuildOscillatedDist(evIt->first, num_dimensions, pdfConfig, dataSet, deltam21_nom, theta12_param, indexDistance, reactorRatio);
        dist.AddPadding();
        fakeDataDist = DistBuilder::BuildOscillatedDist(evIt->first, num_dimensions, pdfConfig, dataSet, deltam21_fd, theta12_fdparam, indexDistance, reactorRatio);
        fakeDataDist.AddPadding();

        // Now we will scale the constraint on the unoscillated reactor flux by the ratio of the oscillated to unoscillated number of events

        double noms_config = noms[evIt->first];
        if (constrMeans.find(evIt->first) != constrMeans.end())
        {
          constrMeans[evIt->first] = constrMeans[evIt->first] * reactorRatio;
          constrSigmas[evIt->first] = constrSigmas[evIt->first] * reactorRatio;
        }
        noms[evIt->first] = noms_config * reactorRatio;
        mins[evIt->first] = mins[evIt->first] * reactorRatio;
        maxs[evIt->first] = maxs[evIt->first] * reactorRatio;
        fdValues[evIt->first] = fdValues[evIt->first] * reactorRatio;

        // Now loop over deltam points and make a new pdf for each
        for (int iDeltaM = 0; iDeltaM < npoints; iDeltaM++)
        {
          double deltam21 = deltam21_min + (double)iDeltaM * (deltam21_max - deltam21_min) / (npoints - 1);
          // Whichever form of theta12 we have, let's make a sin^2(theta12) to hand to the oscillated dist builder
          double theta12_param = theta12_nom;
          if (hasSinTheta12)
            theta12_param = theta12_nom * theta12_nom;
          else if (hasTheta12)
            theta12_param = sin(M_PI * theta12_nom / 180) * sin(M_PI * theta12_nom / 180);
          BinnedED oscDist = DistBuilder::BuildOscillatedDist(evIt->first, num_dimensions, pdfConfig, dataSet, deltam21, theta12_param, indexDistance, reactorRatio);

          oscDist.AddPadding();
          oscDist.Normalise();
          oscPDFsMap[dsIt->first].push_back(oscDist);
          // Apply nominal systematic variables to the oscillated distribution
          for (std::map<std::string, Systematic *>::iterator systIt = systMap[dsIt->first].begin(); systIt != systMap[dsIt->first].end(); ++systIt)
          {
            // If group is "", we apply to all groups
            if (systGroup[systIt->first] == "" || std::find(pdfGroups[dsIt->first].back().begin(), pdfGroups[dsIt->first].back().end(), systGroup[systIt->first]) != pdfGroups[dsIt->first].back().end())
            {
              double norm;
              oscDist = systIt->second->operator()(oscDist, &norm);
            }
          }
          // oscDist is just the reactor PDF for this point in the oscillation scan, so we scale that by the expected reactor rate
          oscDist.Scale(noms_config * reactorRatio);
          // oscDatasets is a vector of oscillated "asimov" distributions, so will be the sum of all scaled pdfs (after applying nominal systematics)
          // We'll Add the rest of the PDFs later
          oscDatasets.at(iDeltaM).Add(oscDist);
        }
        // And do the same for theta
        for (int iTheta = 0; iTheta < npoints; iTheta++)
        {
          double theta12 = theta12_min + (double)iTheta * (theta12_max - theta12_min) / (npoints - 1);
          // Whichever form of theta12 we have, let's make a sin^2(theta12) to hand to the oscillated dist builder
          double theta12_param = theta12;
          if (hasSinTheta12)
            theta12_param = theta12 * theta12;
          else if (hasTheta12)
            theta12_param = sin(M_PI * theta12 / 180) * sin(M_PI * theta12 / 180);
          BinnedED oscDist = DistBuilder::BuildOscillatedDist(evIt->first, num_dimensions, pdfConfig, dataSet, deltam21_nom, theta12_param, indexDistance, reactorRatio);

          oscDist.AddPadding();
          oscDist.Normalise();
          oscPDFsMap[dsIt->first].push_back(oscDist);
          // Apply nominal systematic variables to the oscillated distribution
          for (std::map<std::string, Systematic *>::iterator systIt = systMap[dsIt->first].begin(); systIt != systMap[dsIt->first].end(); ++systIt)
          {
            // If group is "", we apply to all groups
            if (systGroup[systIt->first] == "" || std::find(pdfGroups[dsIt->first].back().begin(), pdfGroups[dsIt->first].back().end(), systGroup[systIt->first]) != pdfGroups[dsIt->first].back().end())
            {
              double norm;
              oscDist = systIt->second->operator()(oscDist, &norm);
            }
          }
          // oscDist is just the reactor PDF for this point in the oscillation scan, so we scale that by the expected reactor rate
          oscDist.Scale(noms_config * reactorRatio);
          // oscDatasets is a vector of oscillated "asimov" distributions, so will be the sum of all scaled pdfs (after applying nominal systematics)
          // We'll Add the rest of the PDFs later
          oscDatasets.at(npoints + iTheta).Add(oscDist);
        }
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
        IO::SaveHistogram(marginalised.GetHistogram(), asimovDistDir + "/" + evIt->first + "_" + dsIt->first + ".root", dist.GetName());
        IO::SaveHistogram(marginalisedPDF.GetHistogram(), pdfDir + "/" + evIt->first + "_" + dsIt->first + ".root", unscaledPDF.GetName());
        // Also scale fake data dist by fake data value
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
        // For the non-reactor PDFs, add to the oscillated asimovs too!
        if (evIt->first.find("reactor_nubar") != std::string::npos)
        {
          for (int iOscPoints = 0; iOscPoints < 2 * npoints; iOscPoints++)
          {
            oscDatasets.at(iOscPoints).Add(dist);
          }
        }
        // Also scale fake data dist by fake data value
        if (isFakeData)
        {
          fakeDataset.Add(fakeDataDist);
          IO::SaveHistogram(fakeDataDist.GetHistogram(), fakedataDistDir + "/" + evIt->first + "_" + dsIt->first + ".root", fakeDataDist.GetName());
        }
      }
    } // End loop over PDFs

    // And save combined histogram (final asimov dataset). This is with everything (including oscillation parameters) nominal. This is what we
    // compare to to calculate the llh at each point in the scans
    IO::SaveHistogram(asimov.GetHistogram(), outDir + "/asimov_" + dsIt->first + ".root", "asimov");
    if (isFakeData)
    {
      IO::SaveHistogram(fakeDataset.GetHistogram(), outDir + "/fakedata_" + dsIt->first + ".root", "fakedata");
    }

    // Now save the datasets for each point in the oscillation scans. These are the same as the above but with an oscillated reactor PDF
    for (int iOscPoints = 0; iOscPoints < npoints; iOscPoints++)
    {
      double deltam21 = deltam21_min + (double)iOscPoints * (deltam21_max - deltam21_min) / (npoints - 1);
      std::stringstream dmOscAsimovName;
      dmOscAsimovName << "oscAsimov_dm21_" << deltam21 << ".root";
      IO::SaveHistogram(oscDatasets.at(iOscPoints).GetHistogram(), outDir + "/" + dmOscAsimovName.str(), "asimov");

      double theta12 = theta12_min + (double)iOscPoints * (theta12_max - theta12_min) / (npoints - 1);
      std::stringstream thOscAsimovName;
      thOscAsimovName << "oscAsimov_th12_" << theta12 << ".root";
      IO::SaveHistogram(oscDatasets.at(iOscPoints + npoints).GetHistogram(), outDir + "/" + thOscAsimovName.str(), "asimov");
    }

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
    dataDists[dsIt->first] = dataDist;

    // Now build the likelihood
    BinnedNLLH &lh = testStats.emplace_back();
    lh.SetBuffer("energy", 8, 20);
    lh.SetBufferAsOverflow(true);
    // Add our data
    lh.SetDataDist(dataDist);
    // Set whether or not to use Beeston Barlow
    lh.SetBarlowBeeston(beestonBarlowFlag);
    // Add the systematics and any prior constraints
    for (std::map<std::string, Systematic *>::iterator systIt = systMap[dsIt->first].begin(); systIt != systMap[dsIt->first].end(); ++systIt)
      lh.AddSystematic(systIt->second, systGroup[systIt->first]);
    // Add our pdfs
    lh.AddPdfs(pdfMap[dsIt->first], pdfGroups[dsIt->first], genRates[dsIt->first], normFittingStatuses[dsIt->first]);

    // Initialise to nominal values
    for (ParameterDict::iterator parIt = mins.begin(); parIt != mins.end(); ++parIt)
    {
      if (std::find(datasetPars[parIt->first].begin(), datasetPars[parIt->first].end(), dsIt->first) == datasetPars[parIt->first].end())
        continue;
      parameterValues[dsIt->first][parIt->first] = noms[parIt->first];
    }

    // And finally bring it all together
    lh.RegisterFitComponents();

    // Set to these initial values
    lh.SetParameters(parameterValues[dsIt->first]);

  } // End loop over datasets

  // Now combine the LLH for each dataset
  ParameterDict allParVals;
  for (std::map<std::string, ParameterDict>::iterator parMapIt = parameterValues.begin(); parMapIt != parameterValues.end(); parMapIt++)
  {
    for (ParameterDict::iterator parIt = parMapIt->second.begin(); parIt != parMapIt->second.end(); parIt++)
    {
      allParVals[parIt->first] = parIt->second;
    }
  }

  std::vector<TestStatistic *> rawllhptrs;
  for (BinnedNLLH &lh : testStats)
  {
    rawllhptrs.push_back(&lh);
  }
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
  {

    fullLLH.SetConstraint(ratioIt->first, constrRatioParName.at(ratioIt->first), ratioIt->second, constrRatioSigmas.at(ratioIt->first));
  }

  fullLLH.RegisterFitComponents();

  // Calculate the nominal LLH
  double nomllh = fullLLH.Evaluate();

  // Setup outfile
  std::string outFileName = outDir + "/llh_scan.root";
  TFile *outFile = new TFile(outFileName.c_str(), "recreate");

  // Loop over (non-oscillation) parameters
  for (ParameterDict::iterator parIt = mins.begin(); parIt != mins.end(); ++parIt)
  {
    // Get param name
    std::string name = parIt->first;

    // Set scan range from this parameter's max and min values
    // We plot x axis as relative to the nominal value, so if this = 0 we set it to 1
    double nom = noms[name];
    if (nom == 0)
      nom = 1;

    // And a bit of jiggery pokery here to guarantee that the nominal value is one of the scan points
    double width = (maxs[name] - mins[name]) / (npoints);
    int numStepsBelowNom = floor((noms[name] - mins[name]) / width);
    int numStepsAboveNom = floor((maxs[name] - noms[name]) / width);

    double min = noms[name] - numStepsBelowNom * width;
    double max = noms[name] + numStepsAboveNom * width;

    // Make histogram for this parameter
    TString htitle = Form("%s, Asimov Rate: %f", labelName[name].c_str(), nom);
    TH1D *hScan = new TH1D((name + "_full").c_str(), (labelName[name] + "_full").c_str(), npoints, (min - (width / 2)) / nom, (max + (width / 2)) / nom);
    hScan->SetTitle(std::string(htitle + ";" + labelName[name] + " (rel. to Asimov); -(ln L_{full})").c_str());

    std::cout << "Scanning for " << name << std::endl;
    // Now loop from min to max in npoint steps
    for (int i = 0; i < npoints; i++)
    {

      if (i % countwidth == 0)
        std::cout << i << "/" << npoints << " (" << double(i) / double(npoints) * 100 << "%)" << std::endl;

      // Set Parameters
      double parval = hScan->GetBinCenter(i + 1) * nom;
      double tempval = allParVals[name];
      allParVals[name] = parval;
      fullLLH.SetParameters(allParVals);

      // Evaluate LLH (later do sample and penalty)
      double llh = fullLLH.Evaluate();
      // Set bin contents
      hScan->SetBinContent(i + 1, llh - nomllh);

      // Return to nominal value
      allParVals[name] = tempval;
    }
    // Write Histos
    hScan->Write();
  }

  // Now we do the same for the oscillation parameters

  // And a bit of jiggery pokery here to guarantee that the nominal value is one of the scan points
  double width = (deltam21_max - deltam21_min) / (npoints);
  int numStepsBelowNom = floor((deltam21_nom - deltam21_min) / width);
  int numStepsAboveNom = floor((deltam21_max - deltam21_nom) / width);
  double min = deltam21_nom - numStepsBelowNom * width;
  double max = deltam21_nom + numStepsAboveNom * width;

  TString htitle = Form("%s, Nom. Value: %f", labelName["deltam21"].c_str(), deltam21_nom);
  TH1D *hDeltam = new TH1D("deltam21_full", (labelName["deltam21"] + "_full").c_str(), npoints, (min - (width / 2)), (max + (width / 2)));
  hDeltam->SetTitle(std::string(htitle + "; " + labelName["deltam21"] + " (eV^{2}); -(ln L_{full})").c_str());

  std::cout << "Scanning for deltam21" << std::endl;
  for (int iDeltaM = 0; iDeltaM < npoints; iDeltaM++)
  {
    // Now build a second likelihood for varying oscillation params
    // If we use the same one we have problems because the most PDFs are shrunk but the reactor one isn't

    std::vector<BinnedNLLH> oscTestStats;
    for (DSMap::iterator dsIt = dsPDFMap.begin(); dsIt != dsPDFMap.end(); ++dsIt)
    {
      BinnedNLLH &osclh = oscTestStats.emplace_back();
      osclh.SetBuffer("energy", 8, 20);
      osclh.SetBufferAsOverflow(true);
      // Add our 'data'
      osclh.SetDataDist(dataDists[dsIt->first]);

      // Set whether or not to use Beeston Barlow
      osclh.SetBarlowBeeston(beestonBarlowFlag);

      // Add the systematics and any prior constraints
      for (std::map<std::string, Systematic *>::iterator systIt = systMap[dsIt->first].begin(); systIt != systMap[dsIt->first].end(); ++systIt)
        osclh.AddSystematic(systIt->second, systGroup[systIt->first]);

      // Get all the PDFs in a new vector with the required oscillated reactor PDF
      std::vector<BinnedED> pdfvec;
      for (int iPDF = 0; iPDF < pdfMap[dsIt->first].size(); iPDF++)
      {
        if (iPDF != reacPDFNum[dsIt->first])
          pdfvec.push_back(pdfMap[dsIt->first].at(iPDF));
        else
          pdfvec.push_back(oscPDFsMap[dsIt->first].at(iDeltaM));
      }
      osclh.AddPdfs(pdfvec, pdfGroups[dsIt->first], genRates[dsIt->first], normFittingStatuses[dsIt->first]);

      if (iDeltaM % countwidth == 0)
        std::cout << iDeltaM << "/" << npoints << " (" << double(iDeltaM) / double(npoints) * 100 << "%)" << std::endl;

      // Bring it all together
      osclh.RegisterFitComponents();
      osclh.SetParameters(parameterValues[dsIt->first]);
    }

    std::vector<TestStatistic *> rawoscllhptrs;
    for (BinnedNLLH &oscllh : oscTestStats)
    {
      rawoscllhptrs.push_back(&oscllh);
    }
    StatisticSum fullOscLLH = Sum(rawoscllhptrs);

    // Add constraints
    std::vector<std::string> corrPairs;
    for (ParameterDict::iterator corrIt = constrCorrs.begin(); corrIt != constrCorrs.end(); ++corrIt)
    {
      fullOscLLH.SetConstraint(corrIt->first, constrMeans.at(corrIt->first), constrSigmas.at(corrIt->first), constrCorrParName.at(corrIt->first),
                               constrMeans.at(constrCorrParName.at(corrIt->first)), constrSigmas.at(constrCorrParName.at(corrIt->first)), corrIt->second);
      corrPairs.push_back(constrCorrParName.at(corrIt->first));
    }
    for (ParameterDict::iterator constrIt = constrMeans.begin(); constrIt != constrMeans.end(); ++constrIt)
    {

      // Only add single parameter constraint if correlation hasn't already been applied
      if (constrCorrs.find(constrIt->first) == constrCorrs.end() && std::find(corrPairs.begin(), corrPairs.end(), constrIt->first) == corrPairs.end())
        fullOscLLH.SetConstraint(constrIt->first, constrIt->second, constrSigmas.at(constrIt->first));
    }
    for (ParameterDict::iterator ratioIt = constrRatioMeans.begin(); ratioIt != constrRatioMeans.end(); ++ratioIt)
    {

      fullOscLLH.SetConstraint(ratioIt->first, constrRatioParName.at(ratioIt->first), ratioIt->second, constrRatioSigmas.at(ratioIt->first));
    }

    fullOscLLH.RegisterFitComponents();

    // Evaluate LLH and set histogram bin content
    double llh = fullOscLLH.Evaluate();
    hDeltam->SetBinContent(iDeltaM + 1, llh - nomllh);
  }
  hDeltam->Write();

  // And a bit of jiggery pokery here to guarantee that the nominal value is one of the scan points
  width = (theta12_max - theta12_min) / (npoints);
  numStepsBelowNom = floor((theta12_nom - theta12_min) / width);
  numStepsAboveNom = floor((theta12_max - theta12_nom) / width);
  min = theta12_nom - numStepsBelowNom * width;
  max = theta12_nom + numStepsAboveNom * width;

  // Repeat for theta
  htitle = Form("%s, Nom. Value: %f", labelName[theta12name].c_str(), theta12_nom);
  TH1D *hTheta12 = new TH1D((theta12name + "_full").c_str(), (labelName[theta12name] + "_nom_full").c_str(), npoints, (min - (width / 2)), (max + (width / 2)));
  hTheta12->SetTitle(std::string(htitle + "; " + labelName[theta12name] + " (^{o}); -(ln L_{full})").c_str());

  std::cout << "Scanning for " << theta12name << std::endl;
  for (int iTheta12 = 0; iTheta12 < npoints; iTheta12++)
  {
    // Now build a second likelihood for varying oscillation params
    // If we use the same one we have problems because the most PDFs are shrunk but the reactor one isn't

    std::vector<BinnedNLLH> oscTestStats;
    for (DSMap::iterator dsIt = dsPDFMap.begin(); dsIt != dsPDFMap.end(); ++dsIt)
    {
      BinnedNLLH &osclh = oscTestStats.emplace_back();
      osclh.SetBuffer("energy", 8, 20);
      osclh.SetBufferAsOverflow(true);
      // Add our 'data'
      osclh.SetDataDist(dataDists[dsIt->first]);

      // Set whether or not to use Beeston Barlow
      osclh.SetBarlowBeeston(beestonBarlowFlag);

      // Add the systematics and any prior constraints
      for (std::map<std::string, Systematic *>::iterator systIt = systMap[dsIt->first].begin(); systIt != systMap[dsIt->first].end(); ++systIt)
        osclh.AddSystematic(systIt->second, systGroup[systIt->first]);

      // Get all the PDFs in a new vector with the required oscillated reactor PDF
      std::vector<BinnedED> pdfvec;
      for (int iPDF = 0; iPDF < pdfMap[dsIt->first].size(); iPDF++)
      {
        if (iPDF != reacPDFNum[dsIt->first])
          pdfvec.push_back(pdfMap[dsIt->first].at(iPDF));
        else
          pdfvec.push_back(oscPDFsMap[dsIt->first].at(iTheta12 + npoints));
      }
      osclh.AddPdfs(pdfvec, pdfGroups[dsIt->first], genRates[dsIt->first], normFittingStatuses[dsIt->first]);

      // And constraints
      std::vector<std::string> corrPairs;
      for (ParameterDict::iterator corrIt = constrCorrs.begin(); corrIt != constrCorrs.end(); ++corrIt)
      {
        // If either parameter doesn't exist for this dataset, move along
        if (std::find(datasetPars[corrIt->first].begin(), datasetPars[corrIt->first].end(), dsIt->first) == datasetPars[corrIt->first].end() ||
            std::find(datasetPars[constrCorrParName.at(corrIt->first)].begin(), datasetPars[constrCorrParName.at(corrIt->first)].end(), dsIt->first) == datasetPars[constrCorrParName.at(corrIt->first)].end())
          continue;
        osclh.SetConstraint(corrIt->first, constrMeans.at(corrIt->first), constrSigmas.at(corrIt->first), constrCorrParName.at(corrIt->first),
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
          osclh.SetConstraint(constrIt->first, constrIt->second, constrSigmas.at(constrIt->first));
      }
      for (ParameterDict::iterator ratioIt = constrRatioMeans.begin(); ratioIt != constrRatioMeans.end(); ++ratioIt)
      {
        // If parameter doesn't exist for this dataset, move along
        if (std::find(datasetPars[ratioIt->first].begin(), datasetPars[ratioIt->first].end(), dsIt->first) == datasetPars[ratioIt->first].end())
          continue;

        osclh.SetConstraint(ratioIt->first, constrRatioParName.at(ratioIt->first), ratioIt->second, constrRatioSigmas.at(ratioIt->first));
      }

      if (iTheta12 % countwidth == 0)
        std::cout << iTheta12 << "/" << npoints << " (" << double(iTheta12) / double(npoints) * 100 << "%)" << std::endl;

      // Bring it all together
      osclh.RegisterFitComponents();
      osclh.SetParameters(parameterValues[dsIt->first]);

    }

    std::vector<TestStatistic *> rawoscllhptrs;
    for (BinnedNLLH &oscllh : oscTestStats)
    {
      rawoscllhptrs.push_back(&oscllh);
    }
    StatisticSum fullOscLLH = Sum(rawoscllhptrs);

    // Add constraints
    std::vector<std::string> corrPairs;
    for (ParameterDict::iterator corrIt = constrCorrs.begin(); corrIt != constrCorrs.end(); ++corrIt)
    {
      fullOscLLH.SetConstraint(corrIt->first, constrMeans.at(corrIt->first), constrSigmas.at(corrIt->first), constrCorrParName.at(corrIt->first),
                               constrMeans.at(constrCorrParName.at(corrIt->first)), constrSigmas.at(constrCorrParName.at(corrIt->first)), corrIt->second);
      corrPairs.push_back(constrCorrParName.at(corrIt->first));
    }
    for (ParameterDict::iterator constrIt = constrMeans.begin(); constrIt != constrMeans.end(); ++constrIt)
    {

      // Only add single parameter constraint if correlation hasn't already been applied
      if (constrCorrs.find(constrIt->first) == constrCorrs.end() && std::find(corrPairs.begin(), corrPairs.end(), constrIt->first) == corrPairs.end())
        fullOscLLH.SetConstraint(constrIt->first, constrIt->second, constrSigmas.at(constrIt->first));

      for (ParameterDict::iterator ratioIt = constrRatioMeans.begin(); ratioIt != constrRatioMeans.end(); ++ratioIt)
      {

        fullOscLLH.SetConstraint(ratioIt->first, constrRatioParName.at(ratioIt->first), ratioIt->second, constrRatioSigmas.at(ratioIt->first));
      }
    }

    fullOscLLH.RegisterFitComponents();

    // Evaluate LLH and set histogram bin content
    double llh = fullOscLLH.Evaluate();
    hTheta12->SetBinContent(iTheta12 + 1, llh - nomllh);
  }
  hTheta12->Write();

  // Close file
  outFile->Close();
  std::cout << "Wrote scan to " << outFile->GetName() << std::endl;

  delete outFile;
}

int main(int argc, char *argv[])
{
  if (argc != 6)
  {
    std::cout << "\nUsage: fixedosc_llhscan <fit_config_file> <eve_config_file> <pdf_config_file> <syst_config_file> <oscgrid_config_file>" << std::endl;
    return 1;
  }

  std::string fitConfigFile(argv[1]);
  std::string eveConfigFile(argv[2]);
  std::string pdfConfigPath(argv[3]);
  std::string systConfigFile(argv[4]);
  std::string oscgridConfigFile(argv[5]);

  fixedosc_llhscan(fitConfigFile, eveConfigFile, pdfConfigPath, systConfigFile, oscgridConfigFile);

  return 0;
}
