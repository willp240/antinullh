#pragma once

#include <map>
#include <memory>
#include <string>
#include <vector>
#include <unordered_map>

#include <Systematic.h>
#include <SystConfig.hh>
#include <DatasetSetup.hh>

namespace antinufit
{

  // Each dataset owns its configured systematics
  using SystsByDataset = std::map<std::string, std::map<std::string, std::unique_ptr<Systematic>>>;

  SystsByDataset BuildSysts(const SystConfig &cfg,
                            const PDFConfig &pdfConfig,
                            const ParameterDict &noms,
                            const DSMap &dsPDFMap,
                            const std::unordered_map<int, double> &indexDistance);

} // namespace antinufit
