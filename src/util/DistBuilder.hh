#ifndef __ANTINUFIT__DistBuilder__
#define __ANTINUFIT__DistBuilder__

// Antinu headers
#include <PDFConfig.hh>
#include <DistFiller.h>
#include <Functions.hh>

// OXO headers
#include <BinnedED.h>
#include <DataSet.h>

// ROOT headers
#include <TRandom3.h>

// C++ headers
#include <unordered_map>

class BinnedED;
class DataSet;
class AxisCollection;

namespace antinufit
{
  class PDFConfig;
  class EventConfig;

  class DistBuilder
  {
  public:
    static BinnedED Build(const std::string &name, const int, const PDFConfig , DataSet *data_);
    static BinnedED Build(const std::string &name, const PDFConfig , DataSet *data_);
    static BinnedED BuildOscillatedDist(const std::string &, const int, const PDFConfig, DataSet*, double, double, std::unordered_map<int, double>, double &);
    static AxisCollection BuildAxes(const PDFConfig &, const int);
    static AxisCollection BuildAxes(const PDFConfig &);
  };
}
#endif
