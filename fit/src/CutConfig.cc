#include <CutConfig.hh>

namespace bbfit{

double
CutConfig::GetValue(){
  return fValue;
}

void 
CutConfig::SetValue(double d_){
  fValue = d_;
}

double
CutConfig::GetValue2(){
  return fValue2;
}

void 
CutConfig::SetValue2(double d_){
  fValue2 = d_;
}

const std::string& 
CutConfig::GetObs() const{
  return fObs;
}
void
CutConfig::SetObs(const std::string& s_){
  fObs = s_;
}

const std::string&
CutConfig::GetType() const{
  return fType;
}

void
CutConfig::SetType(const std::string& t_){
  fType = t_;
}

const std::string&
CutConfig::GetName() const {
  return fName;
}
void
CutConfig::SetName(const std::string& name_){
  fName = name_;
}
  
}
