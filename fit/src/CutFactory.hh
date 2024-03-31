#ifndef __BBFIT__CUTFactory__
#define __BBFIT__CUTFactory__
#include <string>

class Cut;
namespace bbfit{
class CutFactory{
public:
  static Cut* New(const std::string& name, 
		  const std::string& type_, const std::string& obs_, 
		  double value, double value2 = 0); 
  // v2 only needed for box
};
}
#endif
