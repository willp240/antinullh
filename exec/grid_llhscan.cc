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

// ROOT headers
#include <TH1D.h>

// c++ headers
#include <sys/stat.h>

using namespace antinufit;

void grid_llhscan(const std::string &fitConfigFile_,
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

  // Create output directories
  struct stat st = {0};
  if (stat(outDir.c_str(), &st) == -1)
    mkdir(outDir.c_str(), 0700);

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

  // Load up the oscillation probability grids information
  // This is really just to get the reactors json file
  OscGridConfigLoader oscGridLoader(oscGridConfigFile_);
  OscGridConfig oscGridConfig = oscGridLoader.Load();
  std::map<int, OscGrid *> dummyOscGridMap;

  // Read the reactor distance info
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

  // Store the oscillation parameter scan region. We got these from the fit config file, but are about to delete them as they're
  // not used in the OXO llh
  double deltam21_nom = noms["deltam21"];
  double deltam21_min = mins["deltam21"];
  double deltam21_max = maxs["deltam21"];
  double theta12_nom = noms["theta12"];
  double theta12_min = mins["theta12"];
  double theta12_max = maxs["theta12"];

  // Define the number of points
  int npoints = 150;
  int countwidth = double(npoints) / double(5);

  // A parameter could have been defined in the fit config but isn't associated with a pdf or systematic
  // If so ignore that parameter
  ParameterDict::iterator paramIt = noms.begin();
  size_t numParams = noms.size();
  int reacPDFNum;
  for (int iParam = 0; iParam < numParams; iParam++)
  {
    bool iterate = true;
    if (toGet.find(paramIt->first) == toGet.end() && std::find(fullParamNameVec.begin(), fullParamNameVec.end(), paramIt->first) == fullParamNameVec.end())
    {
      std::cout << paramIt->first << " parameter defined in fit config but not in syst or event config. It will be ignored." << std::endl;
      constrSigmas.erase(paramIt->first);
      constrMeans.erase(paramIt->first);
      mins.erase(paramIt->first);
      maxs.erase(paramIt->first);
      paramIt = noms.erase(paramIt);
      if (paramIt != noms.begin())
        paramIt--;
      else
        iterate = false;
    }
    if (iterate)
      paramIt++;
  }

  // Create the individual PDFs and Asimov components
  std::vector<BinnedED> pdfs;
  std::vector<int> genRates;
  std::vector<std::vector<std::string>> pdfGroups;

  // Create the empty full dist
  BinnedED asimov = BinnedED("asimov", systAxes);
  asimov.SetObservables(dataObs);

  // And we're going to create a scaled "asimov" dataset and reactor PDF for each point in the oscillation parameter scans
  // The oscAsimovs aren't actually used in the scan, but we save them as they're useful to look at. But the teststat object
  // should build exact equivalents if everything has gone to plan
  std::vector<BinnedED> oscAsimovs;
  std::vector<BinnedED> oscPDFs;
  for (int iOscPoints = 0; iOscPoints < 2 * npoints; iOscPoints++)
  {
    BinnedED oscAsmv = BinnedED("asimov", systAxes);
    oscAsmv.SetObservables(dataObs);
    oscAsimovs.push_back(oscAsmv);
  }

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

    // If it's the reactor PDF, two things are different from the others. Firstly we use the BuildOscillatedDist rather than just Build.
    // Also, we're going to loop over oscillation scan points to make the reactor's PDF for that point as well
    if (it->first == "reactor_nubar")
    {
      reacPDFNum = pdfs.size();
      // Build the distribution with oscillation parameters at their nominal values
      dist = DistBuilder::BuildOscillatedDist(it->first, num_dimensions, pdfConfig, dataSet, deltam21_nom, theta12_nom, indexDistance);

      // Now loop over deltam points and make a new pdf for each
      for (int iDeltaM = 0; iDeltaM < npoints; iDeltaM++)
      {
        double deltam21 = deltam21_min + (double)iDeltaM * (deltam21_max - deltam21_min) / npoints;
        std::cout << "Loading " << it->second.GetPrunedPath() << " deltam21: " << deltam21 << ", theta12: " << theta12_nom << std::endl;
        BinnedED oscDist = DistBuilder::BuildOscillatedDist(it->first, num_dimensions, pdfConfig, dataSet, deltam21, theta12_nom, indexDistance);
        oscDist.AddPadding(1E-6);
        oscDist.Normalise();
        oscPDFs.push_back(oscDist);
        // Apply nominal systematic variables to the oscillated distribution
        for (std::map<std::string, Systematic *>::iterator systIt = systMap.begin(); systIt != systMap.end(); ++systIt)
        {
          // If group is "", we apply to all groups
          if (systGroup[systIt->first] == "" || std::find(pdfGroups.back().begin(), pdfGroups.back().end(), systGroup[systIt->first]) != pdfGroups.back().end())
          {
            double distInt = oscDist.Integral();
            oscDist = systIt->second->operator()(oscDist);
            oscDist.Scale(distInt);
          }
        }
        // oscDist is just the reactor PDF for this point in the oscillation scan, so we scale that by the expected reactor rate
        oscDist.Scale(noms[it->first]);
        // oscAsimovs is a vector of "asimov" distributions, so will be the sum of all scaled pdfs (after applying nominal systematics)
        // We'll Add the rest of the PDFs later
        oscAsimovs.at(iDeltaM).Add(oscDist);
      }
      // And do the same for theta
      for (int iTheta = 0; iTheta < npoints; iTheta++)
      {
        double theta12 = theta12_min + (double)iTheta * (theta12_max - theta12_min) / npoints;
        std::cout << "Loading " << it->second.GetPrunedPath() << " deltam21: " << deltam21_nom << ", theta12: " << theta12 << std::endl;
        BinnedED oscDist = DistBuilder::BuildOscillatedDist(it->first, num_dimensions, pdfConfig, dataSet, deltam21_nom, theta12, indexDistance);
        oscDist.AddPadding(1E-6);
        oscDist.Normalise();
        oscPDFs.push_back(oscDist);
        // Apply nominal systematic variables to the oscillated distribution
        for (std::map<std::string, Systematic *>::iterator systIt = systMap.begin(); systIt != systMap.end(); ++systIt)
        {
          // If group is "", we apply to all groups
          if (systGroup[systIt->first] == "" || std::find(pdfGroups.back().begin(), pdfGroups.back().end(), systGroup[systIt->first]) != pdfGroups.back().end())
          {
            double distInt = oscDist.Integral();
            oscDist = systIt->second->operator()(oscDist);
            oscDist.Scale(distInt);
          }
        }
        // oscDist is just the reactor PDF for this point in the oscillation scan, so we scale that by the expected reactor rate
        oscDist.Scale(noms[it->first]);
        // oscAsimovs is a vector of "asimov" distributions, so will be the sum of all scaled pdfs (after applying nominal systematics)
        // We'll Add the rest of the PDFs later
        oscAsimovs.at(npoints + iTheta).Add(oscDist);
      }
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
    for (std::map<std::string, Systematic *>::iterator systIt = systMap.begin(); systIt != systMap.end(); ++systIt)
    {

      // If group is "", we apply to all groups
      if (systGroup[systIt->first] == "" || std::find(pdfGroups.back().begin(), pdfGroups.back().end(), systGroup[systIt->first]) != pdfGroups.back().end())
      {
        double distInt = dist.Integral();
        dist = systIt->second->operator()(dist);
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
      // For the non-reactor PDFs, add to the oscillated asimovs too!
      if (it->first != "reactor_nubar")
      {
        for (int iOscPoints = 0; iOscPoints < 2 * npoints; iOscPoints++)
        {
          oscAsimovs.at(iOscPoints).Add(dist);
        }
      }
    }
  } // End loop over PDFs

  // And save combined histogram (final asimov dataset). This is with everything (including oscillation parameters) nominal. This is what we
  // compare to to calculate the llh at each point in the scans
  IO::SaveHistogram(asimov.GetHistogram(), outDir + "/asimov.root", "asimov");
  // Now save the datasets for each point in the oscillation scans. These are the same as the above but with an oscillated reactor PDF
  for (int iOscPoints = 0; iOscPoints < 2 * npoints; iOscPoints++)
  {
    std::stringstream oscasimovname;
    oscasimovname << "oscasimov_" << iOscPoints << ".root";
    IO::SaveHistogram(oscAsimovs.at(iOscPoints).GetHistogram(), outDir + "/" + oscasimovname.str(), "asimov");
  }

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

  // Calculate the nominal LLH
  double nomllh = lh.Evaluate();

  // Setup outfile
  std::string outFileName = outDir + "/llh_scan.root";
  TFile *outFile = new TFile(outFileName.c_str(), "recreate");

  // Loop over (non-oscillation) parameters
  for (ParameterDict::iterator it = mins.begin(); it != mins.end(); ++it)
  {
    // Get param name
    std::string name = it->first;
    std::cout << "Scanning for " << name << std::endl;

    // Set scan range from this parameter's max and min values
    // We plot x axis as relative to the nominal value, so if this = 0 we set it to 1
    double nom = noms[name];
    if (nom == 0)
      nom = 1;
    double min = mins[name];
    double max = maxs[name];
    // Make histogram for this parameter
    TString htitle = Form("%s, Asimov Rate: %f", name.c_str(), nom);
    TH1D *hScan = new TH1D((name + "_full").c_str(), (name + "_full").c_str(), npoints, min / nom, max / nom);
    hScan->SetTitle(std::string(htitle + ";" + name + " (rel. to Asimov); -(ln L_{full})").c_str());

    // Now loop from min to max in npoint steps
    for (int i = 0; i < npoints; i++)
    {

      if (i % countwidth == 0)
        std::cout << i << "/" << npoints << " (" << double(i) / double(npoints) * 100 << "%)" << std::endl;

      // Set Parameters
      double parval = hScan->GetBinCenter(i + 1) * nom;
      double tempval = parameterValues[name];
      parameterValues[name] = parval;
      lh.SetParameters(parameterValues);

      // Evaluate LLH (later do sample and penalty)
      double llh = lh.Evaluate();
      // Set bin contents
      hScan->SetBinContent(i + 1, llh - nomllh);

      // Return to nominal value
      parameterValues[name] = tempval;
    }
    // Write Histos
    hScan->Write();
  }

  // Now we do the same for the oscillation parameters
  TString htitle = Form("%s, Nom. Value: %f", "deltam21", deltam21_nom);
  TH1D *hDeltam = new TH1D("deltam21_full", "deltam21_full", npoints, deltam21_min, deltam21_max);
  hDeltam->SetTitle(std::string(htitle + "; #Delta m^{2}_{21} (eV^{2}); -(ln L_{full})").c_str());
  std::cout << "Scanning for deltam21" << std::endl;
  for (int iDeltaM = 0; iDeltaM < npoints; iDeltaM++)
  {
    // Now build a second likelihood for varying oscillation params
    // If we use the same one we have problems because the most PDFs are shrunk but the reactor one isn't
    BinnedNLLH osclh;
    // Add our 'data'
    osclh.SetDataDist(dataDist);

    // Add the systematics and any prior constraints
    for (std::map<std::string, Systematic *>::iterator it = systMap.begin(); it != systMap.end(); ++it)
      osclh.AddSystematic(it->second, systGroup[it->first]);

    // Get all the PDFs in a new vector with the required oscillated reactor PDF
    std::vector<BinnedED> pdfvec;
    for (int iPDF = 0; iPDF < pdfs.size(); iPDF++)
    {
      if (iPDF != reacPDFNum)
        pdfvec.push_back(pdfs.at(iPDF));
      else
        pdfvec.push_back(oscPDFs.at(iDeltaM));
    }
    osclh.AddPdfs(pdfvec, pdfGroups, genRates);

    // And set any constraints
    for (ParameterDict::iterator it = constrMeans.begin(); it != constrMeans.end(); ++it)
      osclh.SetConstraint(it->first, it->second, constrSigmas.at(it->first));

    if (iDeltaM % countwidth == 0)
      std::cout << iDeltaM << "/" << npoints << " (" << double(iDeltaM) / double(npoints) * 100 << "%)" << std::endl;

    // Bring it all together
    osclh.RegisterFitComponents();
    osclh.SetParameters(parameterValues);

    // Evaluate LLH and set histogram bin content
    double llh = osclh.Evaluate();
    hDeltam->SetBinContent(iDeltaM + 1, llh - nomllh);
  }
  hDeltam->Write();

  // Repeat for theta
  htitle = Form("%s, Nom. Value: %f", "theta12", theta12_nom);
  TH1D *hTheta12 = new TH1D("theta12_nom_full", "theta12_nom_full", npoints, theta12_min, theta12_max);
  hTheta12->SetTitle(std::string(htitle + "; #theta_{12} (^{o}); -(ln L_{full})").c_str());
  std::cout << "Scanning for theta12" << std::endl;
  for (int iTheta12 = 0; iTheta12 < npoints; iTheta12++)
  {
    // Now build a second likelihood for varying oscillation params
    // If we use the same one we have problems because the most PDFs are shrunk but the reactor one isn't
    BinnedNLLH osclh;
    // Add our 'data'
    osclh.SetDataDist(dataDist);

    // Add the systematics and any prior constraints
    for (std::map<std::string, Systematic *>::iterator it = systMap.begin(); it != systMap.end(); ++it)
      osclh.AddSystematic(it->second, systGroup[it->first]);

    // Get all the PDFs in a new vector with the required oscillated reactor PDF
    std::vector<BinnedED> pdfvec;
    for (int iPDF = 0; iPDF < pdfs.size(); iPDF++)
    {
      if (iPDF != reacPDFNum)
        pdfvec.push_back(pdfs.at(iPDF));
      else
        pdfvec.push_back(oscPDFs.at(iTheta12 + npoints));
    }
    osclh.AddPdfs(pdfvec, pdfGroups, genRates);

    // And set any constraints
    for (ParameterDict::iterator it = constrMeans.begin(); it != constrMeans.end(); ++it)
      osclh.SetConstraint(it->first, it->second, constrSigmas.at(it->first));

    if (iTheta12 % countwidth == 0)
      std::cout << iTheta12 << "/" << npoints << " (" << double(iTheta12) / double(npoints) * 100 << "%)" << std::endl;

    // Bring it all together
    osclh.RegisterFitComponents();
    osclh.SetParameters(parameterValues);

    // Evaluate LLH and set histogram bin content
    double llh = osclh.Evaluate();
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
    std::cout << "\nUsage: grid_llhscan <fit_config_file> <eve_config_file> <pdf_config_file> <syst_config_file> <oscgrid_config_file>" << std::endl;
    return 1;
  }

  std::string fitConfigFile(argv[1]);
  std::string eveConfigFile(argv[2]);
  std::string pdfConfigPath(argv[3]);
  std::string systConfigFile(argv[4]);
  std::string oscgridConfigFile(argv[5]);

  grid_llhscan(fitConfigFile, eveConfigFile, pdfConfigPath, systConfigFile, oscgridConfigFile);

  return 0;
}
