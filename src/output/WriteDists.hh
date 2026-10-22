#pragma once

#include <BuildDatasets.hh>
#include <OutputDirs.hh>

#include <IO.h>

namespace antinufit
{
  // Write input PDFs and expected spectra after dataset preparation has finished
  inline void WriteInputDists(const PreparedDatasetMap &datasets,
                              const DatasetSetup &setup,
                              const OutputDirs &dirs)
  {
    for (PreparedDatasetMap::const_iterator datasetIt = datasets.begin();
         datasetIt != datasets.end(); ++datasetIt)
    {
      const std::string &dsName = datasetIt->first;
      const PreparedDataset &dataset = datasetIt->second;
      for (size_t i = 0; i < dataset.pdfs.size(); ++i)
      {
        const std::string pdfName = dataset.pdfs[i].GetName();
        const std::string fileName = pdfName + "_" + dsName + ".root";

        // Normalised component pdf, projected onto the data observables
        BinnedED pdf = ProjectToDataObs(dataset.pdfs[i], dataset.asimov, setup.dataObs);
        IO::SaveHistogram(pdf.GetHistogram(), dirs.pdfDir + "/" + fileName,
                          dataset.pdfs[i].GetName());

        // Nominal expected counts for this component, after systematics and rate scaling
        // ie this event type's component to the Asimov dataset
        IO::SaveHistogram(dataset.nominalComponents.at(i).GetHistogram(),
                          dirs.asimovDistDir + "/" + fileName,
                          dataset.pdfs[i].GetName());

        if (!dataset.fakeComponents.empty())
        {
          // Dist for this component at the fake-data parameter values
          IO::SaveHistogram(dataset.fakeComponents.at(i).GetHistogram(),
                            dirs.fakedataDistDir + "/" + fileName,
                            dataset.pdfs[i].GetName());
        }
      }
      // Sum of the asimov components for this dataset
      IO::SaveHistogram(dataset.asimov.GetHistogram(),
                        dirs.outDir + "/asimov_" + dsName + ".root", "asimov");
      if (!dataset.fakeComponents.empty())
      {
        // Sum of the fake-data expected components for this dataset
        IO::SaveHistogram(dataset.fakeDataset.GetHistogram(),
                          dirs.outDir + "/fakedata_" + dsName + ".root", "fakedata");
      }
    }
  }

  // Write fitted components, their totals, and the selected data distribution
  inline void WritePostfitDists(
      const std::string &dsName,
      const PreparedDataset &dataset,
      const DatasetSetup &setup,
      const std::vector<BinnedED> &components,
      const BinnedED &total,
      const OutputDirs &dirs)
  {
    for (size_t i = 0; i < components.size(); ++i)
    {
      const std::string pdfName = dataset.pdfs.at(i).GetName();
      IO::SaveHistogram(components[i].GetHistogram(),
                        dirs.postfitDistDir + "/" + pdfName + "_" + dsName + ".root");
    }

    IO::SaveHistogram(total.GetHistogram(),
                      dirs.postfitDistDir + "/postfitdist_" + dsName + ".root");
    BinnedED data = ProjectToDataObs(dataset.dataDist, total, setup.dataObs);
    IO::SaveHistogram(data.GetHistogram(),
                      dirs.outDir + "/data_" + dsName + ".root");
  }
}
