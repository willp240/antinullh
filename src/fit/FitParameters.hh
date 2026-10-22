#pragma once

#include <string>
#include <map>
#include <vector>

#include <ParameterDict.h>
#include <FitConfig.hh>
#include <DatasetSetup.hh>
#include <OscParams.hh>
#include <SystConfig.hh>

namespace antinufit
{
  // Parameter values used by the run after pruning and reactor-rate scaling
  struct FitParameters
  {
    ParameterDict mins;
    ParameterDict maxs;
    ParameterDict noms;
    ParameterDict sigmas;
    ParameterDict nbins;
    ParameterDict fdValues;
    std::map<std::string, bool> fixedPars;
    ParameterDict constraintMeans;
    ParameterDict constraintSigmas;
    ParameterDict ratioMeans;
    ParameterDict ratioSigmas;
    std::map<std::string, std::string> ratioParName;
    ParameterDict fracMeans;
    ParameterDict fracSigmas;
    std::map<std::string, std::string> fracParName;
    ParameterDict shapeMeans;
    ParameterDict shapeSigmas;
    std::map<std::string, std::vector<std::string>> shapeParNames;
    std::map<std::string, std::string> shapeFuncName;
    ParameterDict corrs;
    std::map<std::string, std::string> corrParName;
  };

  inline FitParameters MakeFitParameters(const FitConfig &cfg)
  {
    FitParameters p;
    p.mins = cfg.GetMinima();
    p.maxs = cfg.GetMaxima();
    p.noms = cfg.GetNominals();
    p.sigmas = cfg.GetSigmas();
    p.nbins = cfg.GetNBins();
    p.fdValues = cfg.GetFakeDataVals();
    p.fixedPars = cfg.GetFixPars();
    p.constraintMeans = cfg.GetConstrMeans();
    p.constraintSigmas = cfg.GetConstrSigmas();
    p.ratioMeans = cfg.GetConstrRatioMeans();
    p.ratioSigmas = cfg.GetConstrRatioSigmas();
    p.ratioParName = cfg.GetConstrRatioParName();
    p.fracMeans = cfg.GetConstrFracMeans();
    p.fracSigmas = cfg.GetConstrFracSigmas();
    p.fracParName = cfg.GetConstrFracParName();
    p.shapeMeans = cfg.GetConstrShapeMeans();
    p.shapeSigmas = cfg.GetConstrShapeSigmas();
    p.shapeParNames = cfg.GetConstrShapeParNames();
    p.shapeFuncName = cfg.GetConstrShapeFuncName();
    p.corrs = cfg.GetConstrCorrs();
    p.corrParName = cfg.GetConstrCorrParName();
    return p;
  }

  // Initialise parameters to nominal for a dataset
  inline ParameterDict InitialiseDatasetParams(
      const ParameterDict &noms,
      const std::map<std::string, std::vector<std::string>> &datasetPars,
      const std::string &dsName)
  {
    ParameterDict values;
    const std::vector<std::string> &pars = datasetPars.at(dsName);

    for (std::vector<std::string>::const_reference p : pars)
      values[p] = noms.at(p);

    return values;
  }

  using DatasetParameterMap = std::map<std::string, std::vector<std::string>>;

  // Remove unused parameters and constraints, returning the parameters for each dataset
  DatasetParameterMap PruneUnusedParameters(FitParameters &p,
                                            const DatasetSetup &setup,
                                            const SystConfig &systConfig,
                                            const OscParams &osc);

} // namespace antinufit
