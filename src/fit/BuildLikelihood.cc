#include <BuildLikelihood.hh>

namespace antinufit
{
  void ConfigureDatasetLikelihood(BinnedNLLH &lh,
                                  const SystsByDataset &syst,
                                  const SystConfig &systConfig,
                                  const FitConfig &fitConfig,
                                  const PDFConfig &pdfConfig,
                                  const std::string &dsName,
                                  const PreparedDataset &model,
                                  const std::vector<BinnedED> &pdfs,
                                  const std::vector<int> &generatedRates)
  {
      for (const std::string &axis : pdfConfig.GetAxisNames())
      {
        if (!pdfConfig.HasLLHBufferBins(axis))
          continue;
        const std::pair<int, int> bins = pdfConfig.GetLLHBufferBins(axis);
        lh.SetBuffer(axis, bins.first, bins.second);
      }

      lh.SetDataDist(model.dataDist);
      lh.SetBarlowBeeston(fitConfig.GetBeestonBarlow());
      const std::map<std::string, std::unique_ptr<Systematic>> &datasetSystematics =
          syst.at(dsName);
      const std::map<std::string, std::string> groups = systConfig.GetGroup();
      for (std::map<std::string, std::unique_ptr<Systematic>>::const_iterator entry =
               datasetSystematics.begin(); entry != datasetSystematics.end(); ++entry)
        lh.AddSystematic(entry->second.get(), groups.at(entry->first));
      lh.AddPdfs(pdfs, model.pdfGroups, generatedRates,
                 &model.normFittingStatuses);
      lh.RegisterFitComponents();
      lh.SetParameters(model.initialValues);
  }

  std::vector<BinnedNLLH> BuildDatasetLikelihoods(
      const PreparedDatasetMap &datasets,
      const SystsByDataset &syst,
      const SystConfig &systConfig,
      const FitConfig &fitConfig,
      const PDFConfig &pdfConfig)
  {
    std::vector<BinnedNLLH> llhs;
    llhs.reserve(datasets.size());
    for (PreparedDatasetMap::const_iterator entry = datasets.begin();
         entry != datasets.end(); ++entry)
    {
      BinnedNLLH &lh = llhs.emplace_back();
      ConfigureDatasetLikelihood(lh, syst, systConfig, fitConfig, pdfConfig,
                                 entry->first, entry->second,
                                 entry->second.pdfs, entry->second.genRates);
    }
    return llhs;
  }
}
