#include <FitConfigLoader.hh>
#include <OscGridConfigLoader.hh>
#include <DatasetSetup.hh>
#include <OutputDirs.hh>
#include <BuildDatasets.hh>
#include <ParamUtils.hh>
#include <FitInputs.hh>
#include <Rand.h>
#include <cmath>
#include <cctype>
#include <iostream>
#include <iomanip>
#include <map>
#include <vector>

using namespace antinufit;

static bool EndsWith(const std::string &s, const std::string &suffix)
{
  return s.size() >= suffix.size() && s.compare(s.size()-suffix.size(), suffix.size(), suffix) == 0;
}

static std::string PairKey(const std::string &n)
{
  // Classifier parameters use explicit detector suffixes rather than the
  // trailing-2 convention used by the other PPO/BisMSB parameters.
  if (n == "class_a_ppo" || n == "class_a_bismsb") return "class_a";
  if (n == "class_s_ppo" || n == "class_s_bismsb") return "class_s";
  // Fit-config convention: the corresponding BisMSB parameter has a
  // trailing "2", while the PPO parameter does not.
  return EndsWith(n, "2") ? n.substr(0, n.size()-1) : n;
}

static std::string NormalisationName(const std::string &name)
{
  std::string base = EndsWith(name, "2") ? name.substr(0, name.size()-1) : name;
  return base + "_norm";
}

static double ValidScaleValue(const std::string &name, double value, double minimum)
{
  // OXO Scale systematics require a strictly positive factor.  A 100%
  // downward variation of a nominal factor of one would otherwise produce 0.
  if (value <= 0 && name.find("energy_scale") != std::string::npos)
    return minimum > 0 ? minimum : 1.0e-12;
  return value;
}

static double ManualSigma(const std::string &name)
{
  // Temporary analysis-specific overrides for unconstrained parameters.
  if (name == "geonu_U") return 2.90;
  if (name == "geonu_U2") return 5.22;
  if (name == "geonu_Th") return 1.40;
  if (name == "geonu_Th2") return 2.52;
  if (name == "p_recoil_energy_scale") return 0.01;
  if (name == "p_recoil_energy_scale2") return 0.11;
  if (name == "class_a_ppo") return 0.25;
  if (name == "class_a_bismsb") return 0.46;
  return 0.0;
}

int main(int argc, char **argv)
{
  if (argc != 6 && argc != 7) { std::cerr << "Usage: param_syst_scan <fit> <events> <pdf> <syst> <oscgrid>\n"; return 1; }
  Rand::SetSeed(0);
  FitConfig config = FitConfigLoader(argv[1]).LoadActive();
  const FitInputs original = MakeFitInputs(config);
  DatasetSetup setup = LoadDatasetSetup(argv[2], argv[3]);
  const std::unordered_map<int,double> distances = LoadIndexDistanceMap(OscGridConfigLoader(argv[5]).Load().GetReactorsJsonFile());
  FitInputs base = original;
  OscParams osc = BuildOscParams(base.p.noms);
  SystSetup syst = BuildSystematics(argv[4], setup.pdfConfig, base.p.noms, setup.dsPDFMap, distances);
  DatasetParsInfo info = PruneUnusedParameters(base, setup, syst, osc);
  EnsureDir(original.run.outDir);
  EnsureDir(original.run.outDir + "/param_scan");

  // Save an unshifted prefit nominal output alongside the parameter scans.
  // The +/-1 sigma variations below are defined relative to these values.
  {
    FitInputs nominal = base;
    nominal.run.outDir = original.run.outDir + "/param_scan/nominal";
    nominal.run.saveOutputs = true;
    OutputDirs nominalDirs = MakeOutputDirs(nominal.run.outDir);
    EnsureDirs(nominalDirs);
    BuildDatasetsConfig nominalCfg;
    nominalCfg.saveShiftedPdfs = true;
    BuildDatasets(nominal, nominalDirs, setup, syst,
                  BuildOscParams(nominal.p.noms), distances, info, nominalCfg);
  }

  std::map<std::string,std::vector<std::string>> groups;
  for (ParameterDict::const_iterator it=base.p.noms.begin(); it!=base.p.noms.end(); ++it) {
    const std::string &n=it->first;
    // The _norm parameters scale PPO and BisMSB together and are not part
    // of this shape/systematic study.
    if (EndsWith(n, "_norm") || EndsWith(n, "_norm2")) continue;
    groups[PairKey(n)].push_back(n);
  }
  for (std::map<std::string,std::vector<std::string>>::const_iterator g=groups.begin(); g!=groups.end(); ++g)
    for (int sign=-1; sign<=1; sign+=2) {
      std::cout << "\n========================================\n"
                << "Starting " << g->first << (sign>0 ? " +1 sigma" : " -1 sigma") << "\n"
                << "========================================\n";
      FitInputs in=base;
      std::map<std::string, double> shifted;
      std::map<std::string, double> minusValues;
      std::map<std::string, double> plusValues;
      for (const std::string &p : g->second) {
        double sigma=0;
        const double manualSigma = ManualSigma(p);
        const std::string normName = NormalisationName(p);
        const bool isRate = in.p.noms.find(normName) != in.p.noms.end();
        if (manualSigma > 0)
          sigma = manualSigma;
        else if (isRate && in.constraints.sigmas.find(normName) != in.constraints.sigmas.end())
          sigma = std::abs(in.constraints.sigmas[normName]) * std::abs(in.p.noms[p]);
        else if (isRate)
          sigma = std::abs(in.p.noms[p]);
        else if (in.constraints.sigmas.find(p)!=in.constraints.sigmas.end())
          sigma = std::abs(in.constraints.sigmas[p]);
        else if (in.p.noms[p] != 0)
          sigma = std::abs(in.p.noms[p]);
        else
          sigma = sign > 0 ? in.p.maxs[p] - in.p.noms[p] : in.p.noms[p] - in.p.mins[p];

        double value = in.p.noms[p] + sign*sigma;
        if (p.find("conv") != std::string::npos && value < 0) value = 0;
        value = ValidScaleValue(p, value, in.p.mins[p]);
        shifted[p] = value;
        double minusValue = in.p.noms[p] - sigma;
        double plusValue = in.p.noms[p] + sigma;
        if (p.find("conv") != std::string::npos) {
          if (minusValue < 0) minusValue = 0;
          if (plusValue < 0) plusValue = 0;
        }
        minusValue = ValidScaleValue(p, minusValue, in.p.mins[p]);
        plusValue = ValidScaleValue(p, plusValue, in.p.mins[p]);
        minusValues[p] = minusValue;
        plusValues[p] = plusValue;
        in.p.noms[p] = value;
      }
      std::cout << "Parameter values for this variation:\n"
                << std::setprecision(10);
      for (const std::string &p : g->second)
        std::cout << "  " << p << ": nominal " << base.p.noms.at(p)
                  << ", -1 sigma " << minusValues.at(p)
                  << ", +1 sigma " << plusValues.at(p) << "\n";
      in.run.outDir=original.run.outDir+"/param_scan/"+g->first+(sign>0?"_plus1sigma":"_minus1sigma");
      in.run.saveOutputs=true;
      OutputDirs dirs=MakeOutputDirs(in.run.outDir); EnsureDirs(dirs);
      BuildDatasetsConfig cfg; cfg.saveShiftedPdfs=true;
      BuildDatasets(in,dirs,setup,syst,BuildOscParams(in.p.noms),distances,info,cfg);
      std::cout << "\n----------------------------------------\n"
                << "DONE: " << g->first << (sign>0 ? " +1 sigma" : " -1 sigma") << "\n"
                << "Output: " << in.run.outDir
                << "\n----------------------------------------\n" << std::endl;
    }
  return 0;
}
