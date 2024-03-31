#include <string>
#include <FitConfigLoader.hh>
#include <DistConfigLoader.hh>
#include <DistConfig.hh>
#include <DistBuilder.hh>
#include <CutConfigLoader.hh>
#include <FitConfigLoader.hh>
#include <FitConfig.hh>
#include <CutFactory.hh>
#include <CutCollection.h>
#include <SystConfigLoader.hh>
#include <SystConfig.hh>
#include <SystFactory.hh>
#include <fstream>
#include <ROOTNtuple.h>
#include <BinnedNLLH.h>
#include <sys/stat.h>
#include <Rand.h>
#include <AxisCollection.h>
#include <IO.h>
#include <TH1D.h>

#include <Scale.h>

using namespace bbfit;

void
llh_scan(const std::string& mcmcConfigFile_, 
	 const std::string& distConfigFile_,
	 const std::string& cutConfigFile_, 
	 const std::string& systConfigFile_,
	 const std::string& dataPath_,
	 const std::string& asmvRatesPath,
	 const std::string& dims_,
	 const std::string& outDir_){
  Rand::SetSeed(0);


  // Load up the configuration data
  FitConfig mcConfig;
  
  typedef std::vector<CutConfig> CutVec;
  CutVec cutConfs;
    
  {
    FitConfigLoader mcLoader(mcmcConfigFile_);
    mcConfig = mcLoader.LoadActive();
    
    CutConfigLoader cutConfLoader(cutConfigFile_);
    cutConfs = cutConfLoader.LoadActive();
  }

  struct stat st = {0};
  if (stat(outDir_.c_str(), &st) == -1) {
    mkdir(outDir_.c_str(), 0700);
  }
    

  // Make the cuts
  CutCollection cutCol;
  for(CutVec::iterator it = cutConfs.begin(); it != cutConfs.end();
      ++it){
    std::string name = it->GetName();
    std::string type = it->GetType();
    std::string obs = it->GetObs();
    double val = it->GetValue();
    double val2 = it->GetValue2();
    Cut *cut = CutFactory::New(name, type, obs, val, val2);
    cutCol.AddCut(*cut);
    delete cut; // cut col takes its own copy
  }
  

  // Load up the dists
  DistConfigLoader dLoader(distConfigFile_);
  DistConfig pConfig = dLoader.Load();
  std::string distDir = pConfig.GetPDFDir();

  std::vector<BinnedED> dists;
  std::vector<int> genRates;
    
  // the ones you actually want to fit are those listed in mcmcconfig
  typedef std::set<std::string> StringSet;
  StringSet distsToFit = mcConfig.GetParamNames();
  
  for(StringSet::iterator it = distsToFit.begin(); it != distsToFit.end();
      ++it){
    std::string distPath = distDir + "/" + *it + ".h5";
    dists.push_back(BinnedED(*it, IO::LoadHistogram(distPath)));
    std::string rootPath = distDir + "/" + *it + ".root";
    TFile *pdfFile = new TFile(rootPath.c_str(), "READ");
    std::vector<int> nGeneratedEvents;
    std::vector<int>* tempVec;
    pdfFile->GetObject("nGeneratedEvents",tempVec);
    nGeneratedEvents = *tempVec;
    pdfFile->Close();
    genRates.push_back(nGeneratedEvents.at(0));
  }

  // if its a root tree then load it up
  BinnedED dataDist;

  if(dataPath_.substr(dataPath_.find_last_of(".") + 1) == "h5"){
    Histogram loaded = IO::LoadHistogram(dataPath_);
    dataDist = BinnedED("data", loaded);
    dataDist.SetObservables(pConfig.GetBranchNames());
  }
  else{
    // Load up the data set
    ROOTNtuple dataToFit(dataPath_, "pruned");
      
    // Log the effects of the cuts
    CutLog log(cutCol.GetCutNames());
      
    // and bin the data inside
    dataDist = DistBuilder::Build("data", pConfig, (DataSet*)&dataToFit, cutCol, log);
       
    std::ofstream ofs((outDir_ + "/data_cut_log.txt").c_str());
    ofs << "Cut log for data set " << dataPath_ << std::endl;
    ofs << log.AsString() << std::endl;
    ofs.close();
  }

  //marginalise over PSD for 3D fitting
  if(dims_=="3d"){
    std::cout<< "Marginilising for 3d" << std::endl;
    std::vector<std::string> keepObs;
    keepObs.push_back("energy");
    keepObs.push_back("r");
    keepObs.push_back("timePSD");
    dataDist = dataDist.Marginalise(keepObs);
  }

  
  //marginalise over PSD for 2D fitting
  if(dims_=="2d"){
    std::cout<< "Marginilising for 2d" << std::endl;
    std::vector<std::string> keepObs;
    keepObs.push_back("energy");
    keepObs.push_back("r");
    dataDist = dataDist.Marginalise(keepObs);
  }

  // Load up the systematics
  SystConfigLoader systLoader(systConfigFile_);
  SystConfig systConfig = systLoader.LoadActive();

  AxisCollection systAxes = DistBuilder::BuildAxes(pConfig);
  std::vector<std::string> dataObs = pConfig.GetBranchNames();
  ObsSet dataObsSet(dataObs);

  ParameterDict syst_nom =    systConfig.GetNominal();
  ParameterDict syst_maxima = systConfig.GetMaxima();
  ParameterDict syst_minima = systConfig.GetMinima();
  ParameterDict syst_mass =   systConfig.GetMass();
  ParameterDict syst_nbins =  systConfig.GetNBins();
  std::map<std::string, std::string> syst_type = systConfig.GetType();
  std::map<std::string, std::string> syst_obs =  systConfig.GetObs();

  std::vector <Systematic*> syst_vec;
  bool beestonBarlowFlag = mcConfig.GetBeestonBarlow();

  //Loop over systematics and declare each type. Must be a better way to do this but
  // it will do for now
  for(std::map<std::string, std::string>::iterator it = syst_type.begin(); it != syst_type.end()\
	;
      ++it) {
    std::vector<std::string> Obs;
    Obs.push_back(syst_obs[it->first]);
    ObsSet obsSet(Obs);
    
    double stddev_nom = 0;
    if(syst_type[it->first] == "conv")
      stddev_nom = syst_nom[it->first+"_stddevs"];
    Systematic *syst = SystFactory::New(it->first, syst_type[it->first], syst_nom[it->first], stddev_nom);
    syst->SetAxes(systAxes);
    syst->SetTransformationObs(obsSet);
    syst->SetDistributionObs(dataObsSet);
    syst->Construct();
    syst_vec.push_back(syst);
  }
  
  // now build the likelihood
  BinnedNLLH lh;
  lh.SetBufferAsOverflow(true);
  lh.SetBuffer("energy",12,12);
  lh.AddPdfs(dists, genRates);
  lh.SetCuts(cutCol);
  lh.SetDataDist(dataDist);
  lh.SetBarlowBeeston(beestonBarlowFlag);
  for (int i_syst = 0; i_syst<syst_vec.size(); i_syst++) {
    lh.AddSystematic(syst_vec.at(i_syst));
  }

  ParameterDict constrMeans  = mcConfig.GetConstrMeans();
  ParameterDict constrSigmas = mcConfig.GetConstrSigmas();
  ParameterDict mins = mcConfig.GetMinima();
  ParameterDict maxs = mcConfig.GetMaxima();

  for(ParameterDict::iterator it = constrMeans.begin(); it != constrMeans.end();
      ++it)
    lh.SetConstraint(it->first, it->second, constrSigmas.at(it->first));

  //Load in asimov rates for central value of each parameter
  TFile *asmvRatesFile = new TFile(asmvRatesPath.c_str(), "OPEN");
  asmvRatesFile->cd();
  std::map<std::string, double>* tempMap;
  asmvRatesFile->GetObject("AsimovRates",tempMap);
  ParameterDict asimovRates = (ParameterDict)*tempMap;

  for(ParameterDict::iterator it = syst_mass.begin(); it != syst_mass.end(); ++it){
    mins[it->first] = syst_minima[it->first];
    maxs[it->first] = syst_maxima[it->first];
    asimovRates[it->first] = syst_nom[it->first];
  }

  lh.RegisterFitComponents();

  //number of points in scan
  int npoints= 150;
  int countwidth = double(npoints)/double(5);
  
  ParameterDict parameterValues;// = asimovRates;

  for(ParameterDict::iterator it = mins.begin(); it != mins.end(); ++it){
    parameterValues[it->first] = asimovRates[it->first];
  }
  for(ParameterDict::iterator it = constrMeans.begin(); it != constrMeans.end();
      ++it)
    parameterValues[it->first] = constrMeans[it->first];
  lh.SetParameters(parameterValues);

  // Outfile
  std::string outFileName = outDir_ + "/llh_scan.root";
  TFile *outFile = new TFile(outFileName.c_str(),"recreate");

  //Loop over parameters
  for(ParameterDict::iterator it = mins.begin(); it != mins.end(); ++it){
    //Get param name
    std::string name = it->first;
    std::cout << "Scanning for " << name << std::endl;
    double nom = asimovRates[name];
    if(nom==0) nom = 1;
    double min = mins[name];
    double max = maxs[name];
    std::cout << "Scanning for " << name << " from " << min << " to " << max << ". Nom " << nom << std::endl;
    //Make histos
    TString htitle = Form("%s, Asimov Rate: %f", name.c_str(), nom);
    TH1D *hScan = new TH1D((name+"_full").c_str(), (name+"_full").c_str(), npoints, min/nom, max/nom);
    hScan->SetTitle(std::string(htitle + ";" + name + " (rel. to Asimov); -(ln L_{full})").c_str());

    //Commented out bits here are in anticipation of splitting LLH calculation in oxo at some point. Could then see prior and sample contributions separately
    //TH1D *hScanSam = new TH1D((name+"_sam").c_str(), (name+"_sam").c_str(), npoints, min/nom, max/nom);
    //hScanSam->SetTitle(std::string(std::string("2LLH_sam, ") + name + ";" + name + "; -2(ln L_{sample})").c_str());
    //TH1D *hScanPen = new TH1D((name+"_pen").c_str(), (name+"_pen").c_str(), npoints, min/nom, max/nom);
    //hScanPen->SetTitle(std::string(std::string("2LLH_pen, ") + name + ";" + name + "; -2(ln L_{penalty})").c_str());

    //loop from min to max in steps of 150 (might want to do smaller range)
    for(int i=0; i<npoints; i++){

      if (i % countwidth == 0) {
	std::cout << i << "/" << npoints << " (" << double(i)/double(npoints) * 100 << "%)" << std::endl;
      }

      //Set Parameters
      double parval = hScan->GetBinCenter(i+1)*nom;
      double tempval = parameterValues[name];
      parameterValues[name] = parval;
      lh.SetParameters(parameterValues);
      //Eval LLH (later do sample and penalty)
      double LLH = lh.Evaluate();

      //SetBinContents
      hScan->SetBinContent(i+1, LLH);
      //hScanSam->SetBinContent(i+1, LLH);
      //hScanPen->SetBinContent(i+1, LLH);

      //return to asimov value
      parameterValues[name] = tempval;
    }
    //Write Histos
    hScan->Write();
    //hScanSam->Write();
    //hScanPen->Write();
  }
  //Close file
  outFile->Close();
  std::cout << "Wrote scan to " << outFile->GetName() << std::endl;
  
  delete outFile;  
}

int main(int argc, char *argv[]){
  if (argc != 9){
    std::cout << "\nUsage: llh_scan <fit_config_file> <dist_config_file> <cut_config_file> <syst_config_file> <data_to_fit> <asimovRatesFile> <4d,3d or 2d> <outdir>" << std::endl;
    return 1;
  }

  std::string fitConfigFile(argv[1]);
  std::string pdfPath(argv[2]);
  std::string cutConfigFile(argv[3]);
  std::string systConfigFile(argv[4]);
  std::string dataPath(argv[5]);
  std::string asmvRatesPath(argv[6]);
  std::string dims(argv[7]);
  std::string outDir(argv[8]);

  llh_scan(fitConfigFile, pdfPath, cutConfigFile, systConfigFile, dataPath, asmvRatesPath, dims, outDir);

  return 0;
}
