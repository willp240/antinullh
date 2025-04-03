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
#include <IO.h>
#include <Rand.h>
#include <Minuit.h>

// ROOT headers
#include <TH1D.h>
#include <TKey.h>
#include <TMatrixD.h>

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

  std::string scaledDistDir = outDir + "/scaled_dists";
  std::string pdfDir = outDir + "/pdfs";
  if (saveOutputs)
  {
    // Create output directories
    struct stat st = {0};
    if (stat(outDir.c_str(), &st) == -1)
      mkdir(outDir.c_str(), 0700);
    if (stat(scaledDistDir.c_str(), &st) == -1)
      mkdir(scaledDistDir.c_str(), 0700);
    if (stat(pdfDir.c_str(), &st) == -1)
      mkdir(pdfDir.c_str(), 0700);
  }

  // Load up all the event types we want to contribute
  typedef std::map<std::string, EventConfig> EvMap;
  EventConfigLoader loader(evConfigFile_);
  EvMap toGet = loader.LoadActive();

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
  std::vector<std::string> fullParamNameVec;

  // Load up the oscillation probability grids
  OscGridConfigLoader oscGridLoader(oscGridConfigFile_);
  OscGridConfig oscGridConfig = oscGridLoader.Load();
  std::map<int, OscGrid *> dummyOscGridMap;

  // First read the reactor distance info
  std::string reactorjson = oscGridConfig.GetReactorsJsonFile();
  std::unordered_map<int, double> indexDistance = LoadIndexDistanceMap(reactorjson);

  // Loop over systematics and declare each one
  std::map<std::string, Systematic *> systMap;
  for (std::map<std::string, std::string>::iterator it = systType.begin(); it != systType.end(); ++it)
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
    fullParamNameVec.insert(fullParamNameVec.end(), systParamNames[it->first].begin(), systParamNames[it->first].end());

    // Now build the systematic
    Systematic *syst = SystFactory::New(it->first, systType[it->first], systParamNames[it->first], noms, dummyOscGridMap, indexDistance);
    AxisCollection systAxes = DistBuilder::BuildAxes(pdfConfig, systDistObs[it->first].size());
    syst->SetAxes(systAxes);
    // The "dimensions" the systematic applies to
    syst->SetTransformationObs(systTransObs[it->first]);
    // All the "dimensions" of the dataset
    syst->SetDistributionObs(systDistObs[it->first]);
    syst->Construct();
    systMap[it->first] = syst;
  }

  double deltam21 = noms["deltam21"];
  double theta12 = noms["theta12"];

  // A parameter could have been defined in the fit config but isn't associated with a pdf or systematic
  // If so ignore that parameter
  ParameterDict::iterator it = noms.begin();
  size_t numParams = noms.size();
  std::cout << std::endl;
  for (int iParam = 0; iParam < numParams; iParam++)
  {
    bool iterate = true;
    if (toGet.find(it->first) == toGet.end() && std::find(fullParamNameVec.begin(), fullParamNameVec.end(), it->first) == fullParamNameVec.end())
    {
      std::cout << it->first << " parameter defined in fit config but not in syst or event config. It will be ignored." << std::endl;
      constrSigmas.erase(it->first);
      constrMeans.erase(it->first);
      mins.erase(it->first);
      maxs.erase(it->first);
      sigmas.erase(it->first);
      it = noms.erase(it);
      if (it != noms.begin())
        it--;
      else
        iterate = false;
    }
    if (iterate)
      it++;
  }
  std::cout << std::endl;
  PrintParams(noms, mins, maxs, constrMeans, constrSigmas, constrRatioMeans, constrRatioSigmas, constrRatioParName, constrCorrs, constrCorrParName);

  // Create the individual PDFs and Asimov components (could these be maps not vectors?)
  std::vector<BinnedED> pdfs;
  std::vector<int> genRates;
  std::vector<std::vector<std::string>> pdfGroups;
  std::vector<NormFittingStatus> *norm_fitting_statuses = new std::vector<NormFittingStatus>;
  double reactorRatio = 1.0; // Ratio of oscillated to unoscillated number of reactor IBDs

  // Create the empty full dist
  BinnedED asimov = BinnedED("asimov", systAxes);
  asimov.SetObservables(dataObs);
  // And an empty fake data dist
  BinnedED fakeDataset = BinnedED("fake_dataset", systAxes);
  fakeDataset.SetObservables(dataObs);

  // Build each of the PDFs, scale them to the correct size
  for (EvMap::iterator it = toGet.begin(); it != toGet.end(); ++it)
  {

    std::cout << "Building distribution for " << it->first << std::endl;
    pdfGroups.push_back(it->second.GetGroup());

    // Get MC events from file
    DataSet *dataSet;
    try
    {
      std::cout << "Loading " << it->second.GetPrunedPath() << std::endl;
      dataSet = new ROOTNtuple(it->second.GetPrunedPath(), "pruned");
    }
    catch (const IOError &e_)
    {
      std::cout << "Warning: skipping " << it->first << " couldn't open files:\n\t" << e_.what() << std::endl;
      continue;
    }

    // Build distribution of those events
    BinnedED dist;
    BinnedED fakeDataDist;
    int num_dimensions = it->second.GetNumDimensions();

    if (it->first == "reactor_nubar")
    {
      // Build the distribution with oscillation parameters at their nominal values
      dist = DistBuilder::BuildOscillatedDist(it->first, num_dimensions, pdfConfig, dataSet, deltam21, theta12, indexDistance, reactorRatio);

      // Now we will scale the constraint on the unoscillated reactor flux by the ratio of the oscillated to unoscillated number of events
      constrMeans[it->first] = constrMeans[it->first] * reactorRatio;
      constrSigmas[it->first] = constrSigmas[it->first] * reactorRatio;
      noms[it->first] = noms[it->first] * reactorRatio;
      mins[it->first] = mins[it->first] * reactorRatio;
      maxs[it->first] = maxs[it->first] * reactorRatio;
      fakeDataDist = DistBuilder::BuildOscillatedDist(it->first, num_dimensions, pdfConfig, dataSet, fdValues["deltam21"], fdValues["theta12"], indexDistance, reactorRatio);
      fdValues[it->first] = fdValues[it->first] * reactorRatio;
    }
    else
    {
      // For all other PDFs, just use Build
      dist = DistBuilder::Build(it->first, num_dimensions, pdfConfig, dataSet);
      fakeDataDist = DistBuilder::Build(it->first, num_dimensions, pdfConfig, dataSet);
    }

    // Save the generated number of events for Beeston Barlow
    genRates.push_back(dist.Integral());

    // Scale for PDF and add to vector
    if (dist.Integral() && fakeDataDist.Integral())
    {
      dist.Normalise();
      fakeDataDist.Normalise();
    }
    pdfs.push_back(dist);

    norm_fitting_statuses->push_back(INDIRECT);

    // Apply nominal systematic variables
    for (std::map<std::string, Systematic *>::iterator systIt = systMap.begin(); systIt != systMap.end(); ++systIt)
    {
      // If group is "", we apply to all groups
      if (systGroup[systIt->first] == "" || std::find(pdfGroups.back().begin(), pdfGroups.back().end(), systGroup[systIt->first]) != pdfGroups.back().end())
      {
        double distInt = dist.Integral();
        double norm;
        dist = systIt->second->operator()(dist, &norm);

        // Set syst parameter to fake data value, and apply to fake data dist and rescale
        std::set<std::string> systParamNames = systMap[systIt->first]->GetParameterNames();
        for (auto itSystParam = systParamNames.begin(); itSystParam != systParamNames.end(); ++itSystParam)
        {
          systMap[systIt->first]->SetParameter(*itSystParam, fdValues[*itSystParam]);
        }
        systIt->second->Construct();
        fakeDataDist = systIt->second->operator()(fakeDataDist, &norm);
      }
    }

    // Now scale the Asimov component by expected count, and also save pdf as a histo
    dist.Scale(noms[it->first]);
    fakeDataDist.Scale(fdValues[it->first]);
    if (dist.GetNDims() != asimov.GetNDims())
    {
      BinnedED marginalised = dist.Marginalise(dataObs);
      asimov.Add(marginalised);
      if (saveOutputs)
        IO::SaveHistogram(marginalised.GetHistogram(), pdfDir + "/" + it->first + ".root", dist.GetName());
      // Also scale fake data dist by fake data value
      if (isFakeData)
      {
        BinnedED marginalisedFakeData = fakeDataDist.Marginalise(dataObs);
        fakeDataset.Add(marginalisedFakeData);
      }
    }
    else
    {
      asimov.Add(dist);
      if (saveOutputs)
        IO::SaveHistogram(dist.GetHistogram(), pdfDir + "/" + it->first + ".root", dist.GetName());
      // Also scale fake data dist by fake data value
      if (isFakeData)
      {
        fakeDataset.Add(fakeDataDist);
      }
    }
  } // End loop over PDFs
  std::cout << std::endl;

  // And save combined histogram (final asimov dataset). This is with everything (including oscillation parameters) nominal. This is what we
  // compare to to calculate the llh
  if (saveOutputs)
  {
    IO::SaveHistogram(asimov.GetHistogram(), outDir + "/asimov.root", "asimov");
    if (isFakeData)
    {
      IO::SaveHistogram(fakeDataset.GetHistogram(), outDir + "/fakedata.root", "fakedata");
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
    std::string dataPath = fitConfig.GetDatafile();

    // Could be h5 or root file
    if (dataPath.substr(dataPath.find_last_of(".") + 1) == "h5")
    {
      Histogram loaded = IO::LoadHistogram(dataPath);
      dataDist = BinnedED("data", loaded);
      dataDist.SetObservables(pdfConfig.GetBranchNames());
    }
    else // If a root file
    {

      TFile *dataFile = TFile::Open(dataPath.c_str(), "READ");
      if (!dataFile || dataFile->IsZombie())
      {
        std::cerr << "Error opening data file " << dataFile << std::endl;
        throw;
      }

      // Get the list of keys in the file
      TList *keyList = dataFile->GetListOfKeys();
      if (!keyList)
      {
        std::cerr << "Error: No keys found in the datafile " << dataPath << std::endl;
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
          ROOTNtuple dataToFit(dataPath, objectName.c_str());

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

  // Now build the likelihood
  BinnedNLLH lh;
  lh.SetBuffer("energy", 1, 20);
  // Add our data
  lh.SetDataDist(dataDist);
  // Set whether or not to use Beeston Barlow
  lh.SetBarlowBeeston(beestonBarlowFlag);
  // Add the systematics and any prior constraints
  for (std::map<std::string, Systematic *>::iterator it = systMap.begin(); it != systMap.end(); ++it)
    lh.AddSystematic(it->second, systGroup[it->first]);
  // Add our pdfs
  lh.AddPdfs(pdfs, pdfGroups, genRates, norm_fitting_statuses);

  // And constraints
  std::vector<std::string> corrPairs;
  for (ParameterDict::iterator it = constrCorrs.begin(); it != constrCorrs.end(); ++it)
  {
    lh.SetConstraint(it->first, constrMeans.at(it->first), constrSigmas.at(it->first), constrCorrParName.at(it->first),
                     constrMeans.at(constrCorrParName.at(it->first)), constrSigmas.at(constrCorrParName.at(it->first)), it->second);
    corrPairs.push_back(constrCorrParName.at(it->first));
  }
  for (ParameterDict::iterator it = constrMeans.begin(); it != constrMeans.end(); ++it)
  {
    // Only add single parameter constraint if correlation hasn't already been applied
    if (constrCorrs.find(it->first) == constrCorrs.end() && std::find(corrPairs.begin(), corrPairs.end(), it->first) == corrPairs.end())
      lh.SetConstraint(it->first, it->second, constrSigmas.at(it->first));
  }
  for (ParameterDict::iterator it = constrRatioMeans.begin(); it != constrRatioMeans.end(); ++it)
    lh.SetConstraint(it->first, constrRatioParName.at(it->first), it->second, constrRatioSigmas.at(it->first));

  // And finally bring it all together
  lh.RegisterFitComponents();

  // Initialise to nominal values
  ParameterDict parameterValues;
  for (ParameterDict::iterator it = mins.begin(); it != mins.end(); ++it)
    parameterValues[it->first] = noms[it->first];
  // If we have constraints, initialise to those
  for (ParameterDict::iterator it = constrMeans.begin(); it != constrMeans.end(); ++it)
    parameterValues[it->first] = constrMeans[it->first];
  // Set to these initial values
  lh.SetParameters(parameterValues);

  // Now do a fit!
  Minuit min;
  min.SetMethod("Migrad");
  min.SetMaxCalls(10000000);
  min.SetMinima(mins);
  min.SetMaxima(maxs);
  min.SetInitialValues(parameterValues);
  min.SetInitialErrors(sigmas);

  FitResult res = min.Optimise(&lh);
  res.SetPrintPrecision(4);
  res.Print();
  ParameterDict bestFit = res.GetBestFit();
  lh.SetParameters(bestFit);
  bool validFit = res.GetValid();
  double finalLLH = lh.Evaluate();

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
    for (ParameterDict::iterator it = bestFit.begin(); it != bestFit.end(); ++it)
    {
      paramNames.push_back(it->first);
      paramVals.push_back(it->second);
      if (validFit)
      {
        paramErr.push_back(sqrt(covMatrix.GetComponent(paramNames.size() - 1, paramNames.size() - 1)));
        covTMatrixD[paramNames.size() - 1][paramNames.size() - 1] = covMatrix.GetComponent(paramNames.size() - 1, paramNames.size() - 1);
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
    std::cout << "Saved fit result to " << outDir + "/fit_result.txt and " << outDir << "r/fit_result.root" << std::endl;

    // Initialise postfit distributions to same axis as data
    BinnedED postfitDist;
    postfitDist = dataDist;
    postfitDist.Empty();

    // Scale the distributions to the correct heights. They are named the same as their fit parameters
    std::cout << "Saving scaled histograms and data to \n\t" << scaledDistDir << std::endl;

    if (dataDist.GetHistogram().GetNDims() < 3)
    {
      ParameterDict bestFit = res.GetBestFit();
      for (size_t i = 0; i < pdfs.size(); i++)
      {
        std::string name = pdfs.at(i).GetName();
        pdfs[i].Normalise();
        // Apply bestfit systematic variables
        for (std::map<std::string, Systematic *>::iterator systIt = systMap.begin(); systIt != systMap.end(); ++systIt)
        {
          // If group is "", we apply to all groups
          if (systGroup[systIt->first] == "" || std::find(pdfGroups.back().begin(), pdfGroups.back().end(), systGroup[systIt->first]) != pdfGroups.back().end())
          {
            std::set<std::string> systParamNames = systMap[systIt->first]->GetParameterNames();
            for (auto itSystParam = systParamNames.begin(); itSystParam != systParamNames.end(); ++itSystParam)
            {
              systMap[systIt->first]->SetParameter(*itSystParam, bestFit[*itSystParam]);
            }
            double norm;
            systIt->second->Construct();
            pdfs[i] = systIt->second->operator()(pdfs[i], &norm);
          }
        }
        pdfs[i].Scale(bestFit[name]);
        IO::SaveHistogram(pdfs[i].GetHistogram(), scaledDistDir + "/" + name + ".root");
        // Sum all scaled distributions to get full postfit "dataset"
        postfitDist.Add(pdfs[i]);
      }

      IO::SaveHistogram(postfitDist.GetHistogram(), scaledDistDir + "/postfitdist.root");
    }
    else
    {
      ParameterDict bestFit = res.GetBestFit();
      for (size_t i = 0; i < pdfs.size(); i++)
      {
        std::string name = pdfs.at(i).GetName();
        pdfs[i].Normalise();
        // Apply bestfit systematic variables
        for (std::map<std::string, Systematic *>::iterator systIt = systMap.begin(); systIt != systMap.end(); ++systIt)
        {
          // If group is "", we apply to all groups
          if (systGroup[systIt->first] == "" || std::find(pdfGroups.back().begin(), pdfGroups.back().end(), systGroup[systIt->first]) != pdfGroups.back().end())
          {
            std::set<std::string> systParamNames = systMap[systIt->first]->GetParameterNames();
            for (auto itSystParam = systParamNames.begin(); itSystParam != systParamNames.end(); ++itSystParam)
            {
              systMap[systIt->first]->SetParameter(*itSystParam, bestFit[*itSystParam]);
            }
            double norm;
            systIt->second->Construct();
            pdfs[i] = systIt->second->operator()(pdfs[i], &norm);
          }
        }
        pdfs[i].Scale(bestFit[name]);

        std::vector<std::string> keepObs;
        keepObs.push_back("energy");
        pdfs[i] = pdfs[i].Marginalise(keepObs);
        IO::SaveHistogram(pdfs[i].GetHistogram(), scaledDistDir + "/" + name + ".root");
        // Sum all scaled distributions to get full postfit "dataset"
        postfitDist.Add(pdfs[i]);
      }

      IO::SaveHistogram(postfitDist.GetHistogram(), scaledDistDir + "/postfitdist.root");
    }

    // And also save the data
    if (dataDist.GetHistogram().GetNDims() < 3)
    {
      IO::SaveHistogram(dataDist.GetHistogram(), scaledDistDir + "/" + "data.root");
    }
    else
    {
      std::vector<std::string> keepObs;
      keepObs.push_back("energy");
      dataDist = dataDist.Marginalise(keepObs);
      IO::SaveHistogram(dataDist.GetHistogram(), scaledDistDir + "/" + "data.root");
    }
  }
  std::cout << "Fit complete for:" << std::endl;
  std::cout << "deltam: " << deltam21 << std::endl;
  std::cout << "theta: " << theta12 << std::endl;
  std::cout << "LLH: " << finalLLH << std::endl;
  std::cout << "FitValid: " << validFit << std::endl;
  std::cout << "ReactorRatio: " << reactorRatio << std::endl
            << std::endl
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
