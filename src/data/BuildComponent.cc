#include <BuildComponent.hh>

#include <ApplySysts.hh>
#include <BuildDatasets.hh>

namespace antinufit
{
  BinnedED BuildComponent(
      const BinnedED &pdf,
      const std::vector<std::string> &pdfGroups,
      const SystsByDataset &syst,
      const std::map<std::string, std::string> &systGroups,
      const std::string &datasetName,
      const ParameterDict &values,
      double rateScale,
      const BinnedED &dataAxes,
      const std::vector<std::string> &dataObs)
  {
    BinnedED component = ApplySystsByGroups(
        pdf, pdfGroups, syst.at(datasetName), systGroups, values);
    component.Scale(values.at(pdf.GetName()) * rateScale);
    return ProjectToDataObs(component, dataAxes, dataObs);
  }

  std::vector<BinnedED> BuildComponents(
      const std::vector<BinnedED> &pdfs,
      const std::vector<std::vector<std::string>> &pdfGroups,
      const SystsByDataset &syst,
      const SystConfig &systConfig,
      const std::string &datasetName,
      const ParameterDict &values,
      const std::vector<double> &rateScales,
      BinnedED &total,
      const std::vector<std::string> &dataObs)
  {
    std::vector<BinnedED> components;
    const std::map<std::string, std::string> systGroups = systConfig.GetGroup();
    components.reserve(pdfs.size());
    for (size_t i = 0; i < pdfs.size(); ++i)
    {
      BinnedED component = BuildComponent(
          pdfs[i], pdfGroups.at(i), syst, systGroups, datasetName,
          values, rateScales.at(i), total, dataObs);
      total.Add(component);
      components.push_back(component);
    }
    return components;
  }
}
