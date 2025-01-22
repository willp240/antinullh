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
#include <Minuit.h>

// ROOT headers
#include <TH1D.h>

// c++ headers
#include <sys/stat.h>

using namespace antinufit;

void grid_fit(const std::string &fitConfigFile_,
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
  bool beestonBarlowFlag = fitConfig.GetBeestonBarlow();
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

  // Create output directories
  struct stat st = {0};
  if (stat(outDir.c_str(), &st) == -1)
    mkdir(outDir.c_str(), 0700);

  std::string scaledDistDir = outDir + "/scaled_dists";
  if (stat(scaledDistDir.c_str(), &st) == -1)
    mkdir(scaledDistDir.c_str(), 0700);

  // Load up all the event types we want to contribute
  typedef std::map<std::string, EventConfig> EvMap;
  EventConfigLoader loader(evConfigFile_);
  EvMap toGet = loader.LoadActive();

  // Load up the PDF information (skeleton axis details, rather than the distributions themselves)
  PDFConfigLoader pdfLoader(pdfConfigFile_);
  PDFConfig pdfConfig = pdfLoader.Load();
  std::string pdfDir = pdfConfig.GetPDFDir();
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

  // Create the individual PDFs and Asimov components (could these be maps not vectors?)
  std::vector<BinnedED> pdfs;
  std::vector<int> genRates;
  std::vector<std::vector<std::string>> pdfGroups;

  // Create the empty full dist
  BinnedED asimov = BinnedED("asimov", systAxes);
  asimov.SetObservables(dataObs);

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

    if (it->first == "reactor_nubar")
    {
      // Build the distribution with oscillation parameters at their nominal values
      dist = DistBuilder::BuildOscillatedDist(it->first, num_dimensions, pdfConfig, dataSet, deltam21, theta12, indexDistance);
    }
    else
    {
      // For all other PDFs, just use Build
      dist = DistBuilder::Build(it->first, num_dimensions, pdfConfig, dataSet);
    }
    // Add small numbers to avoid 0 probability bins
    dist.AddPadding(1E-6);

    // Save the generated number of events for Beeston Barlow
    genRates.push_back(dist.Integral());

    // Scale for PDF and add to vector
    if (dist.Integral())
      dist.Normalise();
    pdfs.push_back(dist);

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
    if (dist.GetNDims() != asimov.GetNDims())
    {
      BinnedED marginalised = dist.Marginalise(dataObs);
      asimov.Add(marginalised);
      IO::SaveHistogram(marginalised.GetHistogram(), pdfDir + "/" + it->first + ".root", dist.GetName());
    }
    else
    {
      asimov.Add(dist);
      IO::SaveHistogram(dist.GetHistogram(), pdfDir + "/" + it->first + ".root", dist.GetName());
    }
  } // End loop over PDFs

  // And save combined histogram (final asimov dataset). This is with everything (including oscillation parameters) nominal. This is what we
  // compare to to calculate the llh
  IO::SaveHistogram(asimov.GetHistogram(), outDir + "/asimov.root", "asimov");

  // Now let's load up the data
  BinnedED dataDist;
  if (!isAsimov)
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
      dataDist = DistBuilder::Build("data", pdfConfig.GetDataAxisCount(), pdfConfig, (DataSet *)&dataToFit);
    }
  }
  else
    dataDist = asimov;

  // Now build the likelihood
  BinnedNLLH lh;
  lh.SetBuffer("energy",1,14);
  // Add our data
  lh.SetDataDist(dataDist);
  // Set whether or not to use Beeston Barlow
  lh.SetBarlowBeeston(beestonBarlowFlag);
  // Add the systematics and any prior constraints
  for (std::map<std::string, Systematic *>::iterator it = systMap.begin(); it != systMap.end(); ++it)
    lh.AddSystematic(it->second, systGroup[it->first]);
  // Add our pdfs
  lh.AddPdfs(pdfs, pdfGroups, genRates);
  // And constraints
  for (ParameterDict::iterator it = constrMeans.begin(); it != constrMeans.end(); ++it)
    lh.SetConstraint(it->first, it->second, constrSigmas.at(it->first));
  for (ParameterDict::iterator it = constrRatioMeans.begin(); it != constrRatioMeans.end(); ++it)
    lh.SetConstraint(it->first, constrRatioParName.at(it->first),it->second, constrRatioSigmas.at(it->first));

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
  min.SetMethod("Simplex");
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

  // Now save the results
  res.SaveAs(outDir + "/fit_result.txt");
  double finalLLH = lh.Evaluate();
  std::ofstream file(outDir + "/fit_result.txt", std::ios::app);
  file << "\nLLH: " << finalLLH << "\n";
  file.close();
  std::cout << "Saved fit result to " << outDir + "/fit_result.txt" << std::endl;

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
      std::set<std::string> systParamNames = systMap[it->first]->GetParameterNames();
      for (auto itSystParam = systParamNames.begin(); itSystParam != systParamNames.end(); ++itSystParam)
      {
        systMap[it->first]->SetParameter(*itSystParam, bestFit[*itSystParam]);
      }
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
      std::set<std::string> systParamNames = systMap[it->first]->GetParameterNames();
      for (auto itSystParam = systParamNames.begin(); itSystParam != systParamNames.end(); ++itSystParam)
      {
        systMap[it->first]->SetParameter(*itSystParam, bestFit[*itSystParam]);
      }
      postfitDist = systMap[it->first]->operator()(postfitDist);
      postfitDist.Scale(distInt);

      IO::SaveHistogram(postfitDist.GetHistogram(), scaledDistDir + "/postfitdist.root");
    }
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

  std::cout << "Fit complete" << std::endl;
}

int main(int argc, char *argv[])
{
  if (argc != 6)
  {
    std::cout << "\nUsage: grid_fit <fit_config_file> <eve_config_file> <pdf_config_file> <syst_config_file> <oscgrid_config_file>" << std::endl;
    return 1;
  }

  std::string fitConfigFile(argv[1]);
  std::string eveConfigFile(argv[2]);
  std::string pdfConfigPath(argv[3]);
  std::string systConfigFile(argv[4]);
  std::string oscgridConfigFile(argv[5]);

  grid_fit(fitConfigFile, eveConfigFile, pdfConfigPath, systConfigFile, oscgridConfigFile);

  return 0;
}
