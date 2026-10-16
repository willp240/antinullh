#pragma once

#include <map>
#include <string>
#include <unordered_map>
#include <vector>

#include <DistBuilder.hh>
#include <EventConfig.hh>
#include <DatasetSetup.hh>
#include <FitInputs.hh>
#include <OscParams.hh>

namespace antinufit
{

  struct OscillationPoint
  {
    double deltam21 = 0.0;
    double theta12SinSq = 0.0;
  };

  struct BuildComponentConfig
  {
    bool buildScanVariants = false;
    std::vector<OscillationPoint> scanPoints;
  };

  struct ComponentBuildResult
  {
    BinnedED dist;
    BinnedED fakeDist;
    BinnedED unscaledPdf;
    int generated = 0;

    double reactorRatio = 1.0;
    double reactorRatioFD = 1.0;

    // Optional scan-only payload: one oscillated PDF per scan point.
    std::vector<BinnedED> scanDists;
    std::vector<double> scanReactorRatios;
  };

  // Builds dist + fakeDist for one PDF component and applies *no systematics*.
  // Mutates `in.p` for oscillated components (scales nominals/mins/maxs/constraints)
  inline ComponentBuildResult BuildComponentDists(FitInputs &in,
                                                  const DatasetSetup &setup,
                                                  const EventConfig &ev,
                                                  const std::string &pdfName,
                                                  int numDimensions,
                                                  DataSet *dataSet,
                                                  const OscParams &osc,
                                                  const std::unordered_map<int, double> &indexDistance,
                                                  std::map<std::string, double> &reactorRatioOut,
                                                  std::map<std::string, double> &reactorRatioFDOut,
                                                  const BuildComponentConfig &cfg = BuildComponentConfig())
  {
    ComponentBuildResult r;

    if (ev.GetOscillated())
    {
      // Ensure keys exist
      reactorRatioOut[pdfName];
      reactorRatioFDOut[pdfName];

      const double theta12_sinsq = osc.theta12_sinsq();
      const double theta12_fd_sinsq =
          Theta12ToSinSq(osc.theta12name, in.p.fdValues.at(osc.theta12name));

      r.dist = DistBuilder::BuildOscillatedDist(pdfName, numDimensions, setup.pdfConfig,
                                                dataSet, osc.deltam21, theta12_sinsq,
                                                indexDistance, reactorRatioOut[pdfName]);

      r.fakeDist = DistBuilder::BuildOscillatedDist(pdfName, numDimensions, setup.pdfConfig,
                                                    dataSet, in.p.fdValues.at("deltam21"),
                                                    theta12_fd_sinsq,
                                                    indexDistance, reactorRatioFDOut[pdfName]);

      r.reactorRatio = reactorRatioOut[pdfName];
      r.reactorRatioFD = reactorRatioFDOut[pdfName];

      if (cfg.buildScanVariants)
      {
        r.scanDists.reserve(cfg.scanPoints.size());
        r.scanReactorRatios.reserve(cfg.scanPoints.size());

        for (std::vector<OscillationPoint>::const_iterator pIt = cfg.scanPoints.begin();
             pIt != cfg.scanPoints.end(); ++pIt)
        {
          double scanRatio = 1.0;
          BinnedED scanDist = DistBuilder::BuildOscillatedDist(
              pdfName,
              numDimensions,
              setup.pdfConfig,
              dataSet,
              pIt->deltam21,
              pIt->theta12SinSq,
              indexDistance,
              scanRatio);

          scanDist.AddPadding();
          if (scanDist.Integral())
            scanDist.Normalise();

          r.scanDists.push_back(scanDist);
          r.scanReactorRatios.push_back(scanRatio);
        }
      }

      // Scale constraints + parameter bounds exactly like your current code
      if (in.constraints.means.find(pdfName) != in.constraints.means.end())
      {
        in.constraints.means[pdfName] *= r.reactorRatio;
        in.constraints.sigmas[pdfName] *= r.reactorRatio;
      }

      in.p.noms[pdfName] *= r.reactorRatio;
      in.p.mins[pdfName] *= r.reactorRatio;
      in.p.maxs[pdfName] *= r.reactorRatio;

      in.p.fdValues[pdfName] *= r.reactorRatioFD;
    }
    else if (ev.GetFlat())
    {
      r.dist = DistBuilder::BuildFlatDist(pdfName, numDimensions, setup.pdfConfig);
      r.fakeDist = DistBuilder::BuildFlatDist(pdfName, numDimensions, setup.pdfConfig);
    }
    else
    {
      r.dist = DistBuilder::Build(pdfName, numDimensions, setup.pdfConfig, dataSet);
      r.fakeDist = DistBuilder::Build(pdfName, numDimensions, setup.pdfConfig, dataSet);
    }

    r.dist.AddPadding();
    r.fakeDist.AddPadding();

    r.generated = static_cast<int>(r.dist.Integral());

    if (r.dist.Integral() && r.fakeDist.Integral())
    {
      r.dist.Normalise();
      r.fakeDist.Normalise();
    }

    r.unscaledPdf = r.dist; // pre-syst, pre-rate scaling
    return r;
  }

} // namespace antinufit
