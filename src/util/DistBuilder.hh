#ifndef __ANTINUFIT__DistBuilder__
#define __ANTINUFIT__DistBuilder__

// Antinu headers
#include <PDFConfig.hh>
#include <DistFiller.h>

// OXO headers
#include <BinnedED.h>

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
    static AxisCollection BuildAxes(const PDFConfig &, const int);
    static AxisCollection BuildAxes(const PDFConfig &);
  };
}
#endif
