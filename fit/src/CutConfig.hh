#ifndef __BBFIT__CutConfig__
#define __BBFIT__CutConfig__
#include <string>

namespace bbfit{
class CutConfig{
public:
  double GetValue();
  void   SetValue(double);

  double GetValue2();
  void   SetValue2(double);

  const std::string& GetObs() const;
  void SetObs(const std::string&);

  const std::string& GetType() const;
  void SetType(const std::string&);

  const std::string& GetName() const;
  void SetName(const std::string& name_);

private:
  std::string fName;
  std::string fType;
  std::string fObs;
  double fValue;  
  double fValue2;
};
}
#endif
