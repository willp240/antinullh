#ifndef __BBFIT__SystFactory__
#define __BBFIT__SystFactory__
#include <string>

class Systematic;
namespace bbfit{
class SystFactory{
public:
  static Systematic* New(const std::string& name, 
		  const std::string& type_,
		  double value, double value2 = 0); 
 };
}
#endif
