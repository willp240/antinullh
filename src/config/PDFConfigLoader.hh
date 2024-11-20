#ifndef __ANTINUFIT__PDFConfigLoader__
#define __ANTINUFIT__PDFConfigLoader__

// Antinu headers
#include <PDFConfig.hh>

// OXO headers
#include <ConfigLoader.hh>

// c++ headers
#include <algorithm>

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
#endif
