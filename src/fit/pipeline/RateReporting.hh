#pragma once

#include <FitInputs.hh>
#include <PDFConfig.hh>

#include <BinnedED.h>
#include <DistTools.h>
#include <TH1D.h>

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <map>
#include <string>
#include <vector>

namespace antinufit
{
  typedef std::map<std::string, std::map<std::string, double>> RateByDatasetAndPdf;

  inline double IntegrateBufferedBins(const BinnedED &dist,
                                      const PDFConfig &pdfConfig,
                                      const std::vector<std::string> &dataObs)
  {
    std::string bufferAxisName;
    for (std::vector<std::string>::const_iterator obsIt = dataObs.begin();
         obsIt != dataObs.end(); ++obsIt)
    {
      if (pdfConfig.HasLLHBufferBins(*obsIt))
      {
        bufferAxisName = *obsIt;
        break;
      }
    }

    BinnedED projectedDist = dist;
    if (!bufferAxisName.empty() && projectedDist.GetNDims() != 1)
    {
      std::vector<std::string> keepObs;
      keepObs.push_back(bufferAxisName);
      projectedDist = projectedDist.Marginalise(keepObs);
    }

    const Histogram projectedHist = projectedDist.GetHistogram();
    if (projectedHist.GetNDims() != 1)
    {
      return projectedDist.Integral();
    }

    TH1D rootHist = DistTools::ToTH1D(projectedHist);
    if (bufferAxisName.empty())
    {
      return rootHist.Integral();
    }

    const std::pair<int, int> bufferBins = pdfConfig.GetLLHBufferBins(bufferAxisName);
    const int nBins = rootHist.GetNbinsX();
    const int binMin = std::max(1, bufferBins.first);
    const int binMax = std::min(nBins, nBins - bufferBins.second);

    if (binMax < binMin)
    {
      return 0.0;
    }
    return rootHist.Integral(binMin, binMax);
  }

  inline void PrintNominalBestFitSummary(
      const FitInputs &in,
      const ParameterDict &bestFit,
      const RateByDatasetAndPdf &nominalRatesByDataset,
      const RateByDatasetAndPdf &bestFitRatesByDataset)
  {
    struct RateComparison
    {
      double nominal = 0.0;
      double bestFit = 0.0;
    };
    typedef std::map<std::string, RateComparison> RateByDataset;
    typedef std::map<std::string, RateByDataset> RateByParameter;

    RateByParameter ratesByParameter;

    for (RateByDatasetAndPdf::const_iterator dsIt = nominalRatesByDataset.begin();
         dsIt != nominalRatesByDataset.end(); ++dsIt)
    {
      const std::string &dsName = dsIt->first;
      for (std::map<std::string, double>::const_iterator pdfIt = dsIt->second.begin();
           pdfIt != dsIt->second.end(); ++pdfIt)
      {
        ratesByParameter[pdfIt->first][dsName].nominal += pdfIt->second;
      }
    }

    for (RateByDatasetAndPdf::const_iterator dsIt = bestFitRatesByDataset.begin();
         dsIt != bestFitRatesByDataset.end(); ++dsIt)
    {
      const std::string &dsName = dsIt->first;
      for (std::map<std::string, double>::const_iterator pdfIt = dsIt->second.begin();
           pdfIt != dsIt->second.end(); ++pdfIt)
      {
        ratesByParameter[pdfIt->first][dsName].bestFit += pdfIt->second;
      }
    }

    std::cout << "\nParameter Summary (Nominal vs Best Fit)\n";
    std::cout << std::left << std::setw(45) << "Parameter Name"
              << std::right << std::setw(18) << "Nominal Value"
              << std::setw(18) << "Bestfit Value" << std::endl;
    std::cout << std::string(81, '-') << std::endl;
    std::cout << std::defaultfloat << std::setprecision(5) << std::showpoint;

    for (ParameterDict::const_iterator parIt = in.p.noms.begin();
         parIt != in.p.noms.end(); ++parIt)
    {
      const std::string &parName = parIt->first;
      const double nominalValue = parIt->second;
      const double bestFitValue = bestFit.at(parName);

      std::cout << std::left << std::setw(45) << parName
                << std::right << std::setw(18) << nominalValue
                << std::setw(18) << bestFitValue << std::endl;
    }

    std::cout << "\nBest-fit PDF Integrals (Non-buffer region)\n";
    std::cout << std::left << std::setw(45) << "PDF / Dataset"
              << std::right << std::setw(18) << "Asimov Integral"
              << std::right << std::setw(18) << "Bestfit Integral" << std::endl;
    std::cout << std::string(81, '-') << std::endl;

    for (RateByParameter::const_iterator pdfIt = ratesByParameter.begin();
         pdfIt != ratesByParameter.end(); ++pdfIt)
    {
      const std::string &pdfName = pdfIt->first;
      const RateByDataset &datasetRates = pdfIt->second;

      for (RateByDataset::const_iterator rateIt = datasetRates.begin();
           rateIt != datasetRates.end(); ++rateIt)
      {
        const std::string rowName = pdfName + " [" + rateIt->first + "]";
        std::cout << std::left << std::setw(45) << rowName
                  << std::right << std::setw(18) << rateIt->second.nominal
                  << std::right << std::setw(18) << rateIt->second.bestFit << std::endl;
      }
    }

    std::cout << std::endl;
  }
} // namespace antinufit
