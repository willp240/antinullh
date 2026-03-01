#pragma once

#include <stdexcept>
#include <string>

#include <TFile.h>
#include <TKey.h>
#include <TList.h>
#include <TH1.h>
#include <TH1D.h>

#include <ROOTNtuple.h>
#include <DistBuilder.hh>
#include <DistTools.h>

namespace antinufit
{

  inline BinnedED LoadDataDist(const std::string &path,
                               const PDFConfig &pdfConfig,
                               const std::string &distName = "data")
  {
    // Open ROOT file
    TFile f(path.c_str(), "READ");
    if (f.IsZombie())
      throw std::runtime_error("Error opening data file: " + path);

    TList *keyList = f.GetListOfKeys();
    if (!keyList || keyList->GetSize() == 0)
      throw std::runtime_error("No keys found in data file: " + path);

    // Iterate keys and load the first supported object
    TIter nextKey(keyList);
    TKey *key = nullptr;

    while ((key = static_cast<TKey *>(nextKey())))
    {
      const std::string className = key->GetClassName();
      const std::string objectName = key->GetName();

      // TNtuple path
      if (className == "TNtuple")
      {
        ROOTNtuple dataToFit(path, objectName.c_str());

        // Build binned distribution from the ntuple
        return DistBuilder::Build(distName.c_str(),
                                  pdfConfig.GetDataAxisCount(),
                                  pdfConfig,
                                  static_cast<DataSet *>(&dataToFit));
      }

      // Histogram path
      if (className == "TH1D" || className.rfind("TH1", 0) == 0)
      {
        TH1D *h = static_cast<TH1D *>(f.Get(objectName.c_str()));
        if (!h)
          throw std::runtime_error("Failed to load histogram '" + objectName + "' from: " + path);

        Histogram loaded = DistTools::ToHist(*h);
        BinnedED dataDist(distName.c_str(), loaded);

        dataDist.SetObservables(pdfConfig.GetDataBranchNames());
        AxisCollection axes = DistBuilder::BuildAxes(pdfConfig, pdfConfig.GetDataAxisCount());
        dataDist.SetAxes(axes);

        return dataDist;
      }
    }

    throw std::runtime_error("No TNtuple or TH1* found in data file: " + path);
  }

  // Returns nullptr if the file can't be opened.
  // Owns the dataset via unique_ptr so you can't leak it.
  inline std::unique_ptr<DataSet> LoadMCDataset(const std::string &path,
                                                const std::string &treeName = "pruned")
  {
    try
    {
      return std::unique_ptr<DataSet>(new ROOTNtuple(path, treeName.c_str()));
    }
    catch (const IOError &e_)
    {
      std::cout << "Warning: couldn't open MC file:\n\t" << e_.what() << std::endl;
      return nullptr;
    }
  }

} // namespace antinufit
