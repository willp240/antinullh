#include <FitConfigLoader.hh>

namespace antinufit
{

  FitConfigLoader::FitConfigLoader(const std::string &filePath_)
  {
    fPath = filePath_;
  }

  FitConfigLoader::~FitConfigLoader()
  {
    ConfigLoader::Close();
  }

  FitConfig
  FitConfigLoader::LoadActive() const
  {
    FitConfig ret;
    ConfigLoader::Open(fPath);
    int it;
    int burnIn;
    int HMCit;
    int HMCburnIn;
    int nSteps;
    double epsilon;
    double sigmaScale;
    std::string outDir;
    std::string dataSet;
    bool beestonBarlowFlag;
    bool asimovFlag;
    double livetime;
    std::string datafile;

    ConfigLoader::Load("summary", "iterations", it);
    ConfigLoader::Load("summary", "burn_in", burnIn);
    ConfigLoader::Load("summary", "hmc_iterations", HMCit);
    ConfigLoader::Load("summary", "hmc_burn_in", HMCburnIn);
    ConfigLoader::Load("summary", "output_directory", outDir);
    ConfigLoader::Load("summary", "n_steps", nSteps);
    ConfigLoader::Load("summary", "epsilon", epsilon);
    ConfigLoader::Load("summary", "sigma_scale", sigmaScale);
    ConfigLoader::Load("summary", "beeston_barlow", beestonBarlowFlag);
    ConfigLoader::Load("summary", "asimov", asimovFlag);
    ConfigLoader::Load("summary", "datafile", datafile);
    ConfigLoader::Load("summary", "livetime", livetime);

    ret.SetOutDir(outDir);
    ret.SetNSteps(nSteps);
    ret.SetEpsilon(epsilon);
    ret.SetIterations(it);
    ret.SetBurnIn(burnIn);
    ret.SetHMCIterations(HMCit);
    ret.SetHMCBurnIn(HMCburnIn);
    ret.SetSigmaScale(sigmaScale);
    ret.SetBeestonBarlow(beestonBarlowFlag);
    ret.SetAsimov(asimovFlag);
    ret.SetDatafile(datafile);
    ret.SetLivetime(livetime);

    typedef std::set<std::string> StringSet;
    StringSet toLoad;
    ConfigLoader::Load("summary", "fit_dists", toLoad);

    std::string name;
    double nom;
    double min;
    double max;
    double sig;
    double constrMean;
    double constrSigma;
    double constrRatioMean;
    double constrRatioSigma;
    std::string constrRatioParName;
    int nbins;

    if (std::find(toLoad.begin(), toLoad.end(), "all") != toLoad.end())
    {
      toLoad = ConfigLoader::ListSections();
      toLoad.erase("summary");
    }

    for (StringSet::iterator it = toLoad.begin(); it != toLoad.end(); ++it)
    {
      name = *it;
      ConfigLoader::Load(name, "min", min);
      ConfigLoader::Load(name, "max", max);
      ConfigLoader::Load(name, "sig", sig);
      ConfigLoader::Load(name, "nbins", nbins);

      try
      {
        ConfigLoader::Load(name, "nom", nom);
      }
      catch (const std::exception &e)
      {
        nom = max + ((max - min) / 2);
      }

      try
      {
        ConfigLoader::Load(name, "constraint_mean", constrMean);
        ConfigLoader::Load(name, "constraint_sigma", constrSigma);
        ret.AddParameter(name, nom, min, max, sig, nbins, constrMean, constrSigma);
      }
      catch (const ConfigFieldMissing &e_)
      {
        try
        {
          ConfigLoader::Load(name, "constraint_ratiomean", constrRatioMean);
          ConfigLoader::Load(name, "constraint_ratiosigma", constrRatioSigma);
          ConfigLoader::Load(name, "constraint_ratioparname", constrRatioParName);
          ret.AddParameter(name, nom, min, max, sig, nbins,constrRatioMean, constrRatioSigma,constrRatioParName);
        }
        catch(const ConfigFieldMissing &e_)
        {
          ret.AddParameter(name, nom, min, max, sig, nbins);

        }
        
      }
    }

    return ret;
  }

}
