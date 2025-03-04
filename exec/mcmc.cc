// Antinu headers
#include <PDFConfigLoader.hh>
#include <DistBuilder.hh>
#include <FitConfigLoader.hh>
#include <EventConfigLoader.hh>
#include <SystConfigLoader.hh>
#include <SystFactory.hh>
#include <OscGridConfigLoader.hh>

// OXO headers
#include <ROOTNtuple.h>
#include <BinnedNLLH.h>
#include <IO.h>
#include <Rand.h>
#include <MetropolisSampler.h>
#include <MCMC.h>

// ROOT headers
#include <TH1D.h>

// c++ headers
#include <sys/stat.h>

using namespace antinufit;

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
  double sigmaScale = fitConfig.GetSigmaScale();

  // Create output directories
  struct stat st = {0};
  if (stat(outDir.c_str(), &st) == -1)
    mkdir(outDir.c_str(), 0700);

  std::string projDir1D = outDir + "/1dlhproj";
  std::string projDir2D = outDir + "/2dlhproj";
  std::string scaledDistDir = outDir + "/scaled_dists";
  std::string pdfDir = outDir + "/pdfs";
  if (stat(projDir1D.c_str(), &st) == -1)
    mkdir(projDir1D.c_str(), 0700);
  if (stat(projDir2D.c_str(), &st) == -1)
    mkdir(projDir2D.c_str(), 0700);
  if (stat(scaledDistDir.c_str(), &st) == -1)
    mkdir(scaledDistDir.c_str(), 0700);
  if (stat(pdfDir.c_str(), &st) == -1)
    mkdir(pdfDir.c_str(), 0700);

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
  std::string outfilename = oscGridConfig.GetFilename();
  std::string reactorjson = oscGridConfig.GetReactorsJsonFile();
  std::map<int, OscGrid *> oscGridMap;

  // First read the reactor distance info
  std::unordered_map<int, double> indexDistance = LoadIndexDistanceMap(reactorjson);
  for (std::unordered_map<int, double>::iterator it = indexDistance.begin(); it != indexDistance.end(); ++it)
  {
    std::string oscGridFileName = outfilename + "_" + std::to_string(it->first) + ".root";
    OscGrid *oscGrid = new OscGrid(oscGridFileName);
    oscGridMap[it->first] = oscGrid;
  }

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
    Systematic *syst = SystFactory::New(it->first, systType[it->first], systParamNames[it->first], noms, oscGridMap, indexDistance);
    AxisCollection systAxes = DistBuilder::BuildAxes(pdfConfig, systDistObs[it->first].size());
    syst->SetAxes(systAxes);
    // The "dimensions" the systematic applies to
    syst->SetTransformationObs(systTransObs[it->first]);
    // All the "dimensions" of the dataset
    syst->SetDistributionObs(systDistObs[it->first]);
    syst->Construct();
    systMap[it->first] = syst;
  }

  // A parameter could have been defined in the fit config but isn't associated with a pdf or systematic
  // If so ignore that parameter
  ParameterDict::iterator it = noms.begin();
  size_t numParams = noms.size();
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
      nbins.erase(it->first);
      it = noms.erase(it);
      if (it != noms.begin())
        it--;
      else
        iterate = false;
    }
    if (iterate)
      it++;
  }

  PrintParams(noms, mins, maxs, constrMeans, constrSigmas, constrRatioMeans, constrRatioSigmas, constrRatioParName, constrCorrs, constrCorrParName);

  // Create the individual PDFs and Asimov components (could these be maps not vectors?)
  std::vector<BinnedED> pdfs;
  std::vector<int> genRates;
  std::vector<std::vector<std::string>> pdfGroups;
  std::vector<NormFittingStatus> *norm_fitting_statuses = new std::vector<NormFittingStatus>;

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
    int num_dimensions = it->second.GetNumDimensions();
    dist = DistBuilder::Build(it->first, num_dimensions, pdfConfig, dataSet);

    // Save the generated number of events for Beeston Barlow
    genRates.push_back(dist.Integral());

    // Scale for PDF and add to vector
    if (dist.Integral())
      dist.Normalise();
    pdfs.push_back(dist);
    norm_fitting_statuses->push_back(INDIRECT);

    // Now make a fake data dist for the event type
    BinnedED fakeDataDist = dist;

    // Apply nominal systematic variables
    for (std::map<std::string, Systematic *>::iterator it = systMap.begin(); it != systMap.end(); ++it)
    {
      // If group is "", we apply to all groups
      if (systGroup[it->first] == "" || std::find(pdfGroups.back().begin(), pdfGroups.back().end(), systGroup[it->first]) != pdfGroups.back().end())
      {
        double distInt = dist.Integral();
        dist = it->second->operator()(dist);
        dist.Scale(distInt);
      }
    }

    // Now scale the Asimov component by expected count, and also save pdf as a histo
    dist.Scale(noms[it->first]);
    fakeDataDist.Scale(fdValues[it->first]);
    if (dist.GetNDims() != asimov.GetNDims())
    {
      BinnedED marginalised = dist.Marginalise(dataObs);
      asimov.Add(marginalised);
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
      IO::SaveHistogram(dist.GetHistogram(), pdfDir + "/" + it->first + ".root", dist.GetName());
      // Also scale fake data dist by fake data value
      if (isFakeData)
      {
        fakeDataset.Add(fakeDataDist);
      }
    }
  }

  // And save combined histogram (final asimov dataset). This is with everything (including oscillation parameters) nominal. This is what we
  // compare to to calculate the llh
  IO::SaveHistogram(asimov.GetHistogram(), outDir + "/asimov.root", "asimov");
  if (isFakeData)
  {
    IO::SaveHistogram(fakeDataset.GetHistogram(), outDir + "/fakedata.root", "fakedata");
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
    std::string dataPath = fitConfig.GetDatafile();

    // Could be h5 or root file
    if (dataPath.substr(dataPath.find_last_of(".") + 1) == "h5")
    {
      Histogram loaded = IO::LoadHistogram(dataPath);
      dataDist = BinnedED("data", loaded);
      dataDist.SetObservables(pdfConfig.GetBranchNames());
    }
    else
    {
      // Load up the data set
      ROOTNtuple dataToFit(dataPath, "pruned");

      // And bin the data inside
      dataDist = DistBuilder::Build("data", pdfConfig, (DataSet *)&dataToFit);
    }
  }

  // Now build the likelihood
  BinnedNLLH lh;
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
    if (!constrCorrs[it->first] && std::find(corrPairs.begin(), corrPairs.end(), it->first) == corrPairs.end())
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

  // And now the optimiser
  // Create something to do sampling
  MetropolisSampler sampler;

  // Multiply parameter sigmas by global sigma scales
  for (ParameterDict::iterator it = sigmas.begin(); it != sigmas.end(); ++it)
    sigmas[it->first] *= sigmaScale;

  sampler.SetSigmas(sigmas);

  MCMC mh(sampler);

  bool saveChain = true;
  mh.SetSaveChain(saveChain);
  mh.SetMaxIter(nsteps);
  mh.SetBurnIn(burnin);
  mh.SetMinima(mins);
  mh.SetMaxima(maxs);

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

  // Lez Do Dis
  const FitResult &res = mh.Optimise(&lh);
  lh.SetParameters(res.GetBestFit());

  // Now save the results
  res.SaveAs(outDir + "/fit_result.txt");
  std::cout << "Saved fit result to " << outDir + "/fit_result.txt" << std::endl;

  MCMCSamples samples = mh.GetSamples();

  if (saveChain)
  {
    // Messy way to make output filename. OutDir is full path to output result directory.
    // Using find_last_of / and stripping everything before it to get directory name (tempString2).
    // Similarly find name of directory above (where pdfs and fake data is saved). Output filename is then aboveDirectory_resultDirectory.root, inside OutDir.
    // Could surely do this in fewer lines, or just call it outputTree.root or something, but nice to have a more unique name
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
  }

  // Save the histograms
  typedef std::map<std::string, Histogram> HistMap;
  const HistMap &proj1D = samples.Get1DProjections();
  const HistMap &proj2D = samples.Get2DProjections();

  std::cout << "Saving LH projections to \n\t" << projDir1D << "\n\t" << projDir2D << std::endl;

  for (HistMap::const_iterator it = proj1D.begin(); it != proj1D.end(); ++it)
    IO::SaveHistogram(it->second, projDir1D + "/" + it->first + ".root");

  for (HistMap::const_iterator it = proj2D.begin(); it != proj2D.end(); ++it)
    IO::SaveHistogram(it->second, projDir2D + "/" + it->first + ".root");

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
      pdfs[i].Scale(bestFit[name]);
      IO::SaveHistogram(pdfs[i].GetHistogram(), scaledDistDir + "/" + name + ".root");
      // Sum all scaled distributions to get full postfit "dataset"
      postfitDist.Add(pdfs[i]);
    }

    for (std::map<std::string, Systematic *>::iterator it = systMap.begin(); it != systMap.end(); ++it)
    {
      double distInt = postfitDist.Integral();
      systMap[it->first]->SetParameter(it->first, bestFit[it->first]);
      postfitDist = systMap[it->first]->operator()(postfitDist);
      postfitDist.Scale(distInt);
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
      pdfs[i].Scale(bestFit[name]);

      std::vector<std::string> keepObs;
      keepObs.push_back("energy");
      pdfs[i] = pdfs[i].Marginalise(keepObs);
      IO::SaveHistogram(pdfs[i].GetHistogram(), scaledDistDir + "/" + name + ".root");
      // Sum all scaled distributions to get full postfit "dataset"
      postfitDist.Add(pdfs[i]);
    }

    for (std::map<std::string, Systematic *>::iterator it = systMap.begin(); it != systMap.end(); ++it)
    {
      double distInt = postfitDist.Integral();
      systMap[it->first]->SetParameter(it->first, bestFit[it->first]);
      postfitDist = systMap[it->first]->operator()(postfitDist);
      postfitDist.Scale(distInt);
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
    std::cout << "\nUsage: llh_scan <fit_config_file> <eve_config_file> <pdf_config_file> <syst_config_file> <oscgrid_config_file>" << std::endl;
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
