#include <DistBuilder.hh>
#include <DistConfig.hh>
#include <BinnedED.h>
#include <DistFiller.h>
#include <DataSet.h>


namespace bbfit{

AxisCollection
DistBuilder::BuildAxes(const DistConfig& config_){
  // Build the axes
  AxisCollection axes;
  
  double min;
  double max;
  int    nBins;
  std::string name;
  std::string texName;
  std::string branchName;

  // save this last one for later
  for(int i = 0; i < config_.GetAxisCount(); i++){
    config_.GetAxis(i, name, branchName, texName, nBins, min, max);
    axes.AddAxis(BinAxis(name, min, max, nBins, texName));
  }
  
  return axes;
}

BinnedED
DistBuilder::Build(const std::string& name_, const DistConfig& pdfConfig_, DataSet* data_, const CutCollection& cuts_, CutLog& log_){
  // Create the axes
  AxisCollection axes = BuildAxes(pdfConfig_);
  
  // put it all together and return
  BinnedED dist(name_, axes);
  dist.SetObservables(pdfConfig_.GetBranchNames());

  // fill it up
  DistFiller::FillDist(dist, *data_, cuts_, log_);

  return dist;
}

}

