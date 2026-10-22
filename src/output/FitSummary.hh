#pragma once

#include <map>
#include <string>
#include <vector>

#include <FitParameters.hh>
#include <OscParams.hh>
#include <PDFConfig.hh>

namespace antinufit
{

  // Combine vector of strings to one string
  inline std::string Join(const std::vector<std::string> &v, const std::string &sep = ", ")
  {
    std::string out;
    for (size_t i = 0; i < v.size(); ++i)
    {
      out += v[i];
      if (i + 1 < v.size())
        out += sep;
    }
    return out;
  }

  void PrintFitSummary(const FitParameters &p,
                       const std::map<std::string, std::vector<std::string>> &datasets,
                       const OscParams &osc,
                       const PDFConfig &pdfConfig);

} // namespace antinufit
