#ifndef __ANTINUFIT__DistBuilder__
#define __ANTINUFIT__DistBuilder__

// Antinu headers
#include <DistConfig.hh>
#include <DistFiller.h>

// OXO headers
#include <BinnedED.h>

class BinnedED;
class DataSet;
class AxisCollection;

namespace antinufit
{
  class DistConfig;
  class EventConfig;

  class DistBuilder
  {
  public:
    static BinnedED Build(const std::string &name, const int, const DistConfig , DataSet *data_);
    static BinnedED Build(const std::string &name, const DistConfig , DataSet *data_);
    static AxisCollection BuildAxes(const DistConfig &, const int);
    static AxisCollection BuildAxes(const DistConfig &);
  };
}
#endif
