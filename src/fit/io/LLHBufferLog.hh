#pragma once

#include <PDFConfig.hh>

#include <iostream>
#include <string>
#include <utility>
#include <vector>

namespace antinufit
{
  inline void PrintLLHBufferConfig(const PDFConfig &pdfConfig)
  {
    std::cout << "LLH buffer bins by axis:" << std::endl;
    const std::vector<std::string> &axisNames = pdfConfig.GetAxisNames();
    for (std::vector<std::string>::const_iterator axisIt = axisNames.begin();
         axisIt != axisNames.end(); ++axisIt)
    {
      const std::string &axisName = *axisIt;
      if (pdfConfig.HasLLHBufferBins(axisName))
      {
        const std::pair<int, int> bins = pdfConfig.GetLLHBufferBins(axisName);
        std::cout << "  " << axisName << ": low=" << bins.first
                  << ", high=" << bins.second << std::endl;
      }
      else
      {
        std::cout << "  " << axisName << ": none" << std::endl;
      }
    }
  }
}
