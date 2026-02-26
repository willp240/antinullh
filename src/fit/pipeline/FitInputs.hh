#pragma once

#include <FitConstraints.hh>

#include <FitConfig.hh>

namespace antinufit
{

  struct FitRunSettings
  {
    bool isAsimov = false;
    bool isFakeData = false;
    bool beestonBarlowFlag = false;
    bool saveOutputs = false;

    double tolerance = 0.0;
    int strategy = 1;
    std::string method;
    int iterations = 0;

    std::string outDir;
  };

  struct FitParamState
  {
    ParameterDict mins;
    ParameterDict maxs;
    ParameterDict noms;
    ParameterDict sigmas;

    ParameterDict fdValues;
    std::map<std::string, bool> fixedPars;
  };

  struct FitInputs
  {
    FitRunSettings run;
    FitParamState p;
    FitConstraints constraints;
  };

  inline FitInputs MakeFitInputs(const FitConfig &cfg)
  {
    FitInputs in;

    // Run settings
    in.run.isAsimov = cfg.GetAsimov();
    in.run.isFakeData = cfg.GetFakeData();
    in.run.beestonBarlowFlag = cfg.GetBeestonBarlow();
    in.run.saveOutputs = cfg.GetSaveOutputs();
    in.run.tolerance = cfg.GetMinuitTolerance();
    in.run.strategy = cfg.GetMinuitStrategy();
    in.run.method = cfg.GetMinuitMethod();
    in.run.iterations = cfg.GetIterations();
    in.run.outDir = cfg.GetOutDir();

    // Parameter state
    in.p.mins = cfg.GetMinima();
    in.p.maxs = cfg.GetMaxima();
    in.p.noms = cfg.GetNominals();
    in.p.sigmas = cfg.GetSigmas();

    in.p.fdValues = cfg.GetFakeDataVals();
    in.p.fixedPars = cfg.GetFixPars();

    // Constraints bundle
    in.constraints.means = cfg.GetConstrMeans();
    in.constraints.sigmas = cfg.GetConstrSigmas();

    in.constraints.ratioMeans = cfg.GetConstrRatioMeans();
    in.constraints.ratioSigmas = cfg.GetConstrRatioSigmas();
    in.constraints.ratioParName = cfg.GetConstrRatioParName();

    in.constraints.fracMeans = cfg.GetConstrFracMeans();
    in.constraints.fracSigmas = cfg.GetConstrFracSigmas();
    in.constraints.fracParName = cfg.GetConstrFracParName();

    in.constraints.shapeMeans = cfg.GetConstrShapeMeans();
    in.constraints.shapeSigmas = cfg.GetConstrShapeSigmas();
    in.constraints.shapeParNames = cfg.GetConstrShapeParNames();
    in.constraints.shapeFuncName = cfg.GetConstrShapeFuncName();

    in.constraints.corrs = cfg.GetConstrCorrs();
    in.constraints.corrParName = cfg.GetConstrCorrParName();

    return in;
  }

} // namespace antinufit
