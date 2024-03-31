#include <SystFactory.hh>
#include <Shift.h>
#include <Scale.h>
#include <Convolution.h>
#include <Gaussian.h>
#include <Exceptions.h>

namespace bbfit{
Systematic*
SystFactory::New(const std::string& name_, 
		const std::string& type_,
		double value, double value2){
  Systematic* syst;

  if(type_ == "scale"){
    Scale* scale = new Scale("scale");
    scale->RenameParameter("scaleFactor",name_);
    scale->SetScaleFactor( value );
    syst = scale;
  }

  else if(type_ == "shift"){
    Shift* shift = new Shift("shift");
    shift->RenameParameter("shift",name_);
    shift->SetShift( value );
    syst = shift;
  }

  else if(type_ == "conv"){
    Convolution* conv = new Convolution("conv");
    Gaussian* gaus = new Gaussian(value, value2, name_);
    gaus->RenameParameter("means_0", name_);
    gaus->RenameParameter("stddevs_0", name_+"_stddevs");
    conv->SetFunction(gaus);
    syst = conv;
  }

  else{
    throw ValueError("Unknown systematic type: " + type_);
  }

  return syst;
}

}
