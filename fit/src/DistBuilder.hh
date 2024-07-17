#ifndef __ANTINUFIT__DistBuilder__
#define __ANTINUFIT__DistBuilder__
#include <string>

class BinnedED;
class DataSet;
class AxisCollection;

namespace antinufit{
class DistConfig;
class EventConfig;

class DistBuilder{
public:
    static BinnedED Build(const std::string& name, const DistConfig&, DataSet* data_);
  static AxisCollection BuildAxes(const DistConfig&);

};
}
#endif
