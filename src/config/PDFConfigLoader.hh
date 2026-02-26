#pragma once

// Antinu headers
#include <PDFConfig.hh>

// OXO headers
#include <ConfigLoader.hh>

namespace antinufit
{
  class PDFConfig;
  class PDFConfigLoader
  {
  public:
    PDFConfigLoader(const std::string &filePath_);
    ~PDFConfigLoader();
    PDFConfig Load() const;

  private:
    std::string fPath;
  };
}
