#include <CutFactory.hh>
#include <BoolCut.h>
#include <BoxCut.h>
#include <LineCut.h>
#include <Exceptions.h>

namespace bbfit{
Cut*
CutFactory::New(const std::string& name_, 
		const std::string& type_, const std::string& obs_,
		double value, double value2){
  Cut* cut;
  if(type_ == "bool" || type_ == "=="){
    cut = new BoolCut(name_, obs_, value);
  }

  else if(type_ == "box"){
    cut = new BoxCut(name_, obs_, value, value2);
  }

  else if (type_ == "line"){
    std::string sidedness;
    if (value2 > 0)
      sidedness = "lower";
    else
      sidedness = "upper";

    cut = new LineCut(name_, obs_, value, sidedness);
  }
  else{
    throw ValueError("Unknown cut type: " + type_);
  }

  return cut;
}

}
