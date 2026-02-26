#pragma once

#include <string>
#include <sys/stat.h>
#include <sys/types.h>

namespace antinufit
{

  struct OutputDirs
  {
    std::string outDir;
    std::string pdfDir;
    std::string asimovDistDir;
    std::string fakedataDistDir;
    std::string postfitDistDir;
  };

  inline OutputDirs MakeOutputDirs(const std::string &outDir)
  {
    OutputDirs d;
    d.outDir = outDir;
    d.pdfDir = outDir + "/unscaled_pdfs";
    d.asimovDistDir = outDir + "/asimov_dists";
    d.fakedataDistDir = outDir + "/fakedata_dists";
    d.postfitDistDir = outDir + "/postfit_dists";
    return d;
  }

  inline void EnsureDir(const std::string &path, mode_t mode = 0700)
  {
    struct stat st = {0};
    if (stat(path.c_str(), &st) == -1)
      mkdir(path.c_str(), mode);
  }

  inline void EnsureDirs(const OutputDirs &d)
  {
    EnsureDir(d.outDir);
    EnsureDir(d.pdfDir);
    EnsureDir(d.asimovDistDir);
    EnsureDir(d.fakedataDistDir);
    EnsureDir(d.postfitDistDir);
  }

} // namespace antinufit
