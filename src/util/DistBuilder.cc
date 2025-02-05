#include <DistBuilder.hh>
#include <iostream>

namespace antinufit
{

  AxisCollection
  DistBuilder::BuildAxes(const PDFConfig &config_)
  {

    int numAxes = config_.GetAxisCount();
    return BuildAxes(config_, numAxes);
  }

  AxisCollection
  DistBuilder::BuildAxes(const PDFConfig &config_, const int numDimensions)
  {

    if (numDimensions > config_.GetAxisCount())
    {
      throw ValueError("Number of axes required is fewer than number declared in config");
    }

    // Build the axes
    AxisCollection axes;

    double min;
    double max;
    int nBins;
    std::string name;
    std::string texName;
    std::string branchName;

    // save this last one for later
    for (int i = 0; i < numDimensions; i++)
    {
      config_.GetAxis(i, name, branchName, texName, nBins, min, max);
      axes.AddAxis(BinAxis(name, min, max, nBins, texName));
    }

    return axes;
  }

  BinnedED
  DistBuilder::Build(const std::string &name_, const PDFConfig pdfConfig_, DataSet *data_)
  {

    int numAxes = pdfConfig_.GetAxisCount();
    return Build(name_, numAxes, pdfConfig_, data_);
  }

  BinnedED
  DistBuilder::Build(const std::string &name_, const int numDimensions_, const PDFConfig pdfConfig_, DataSet *data_)
  {
    // Create the axes
    AxisCollection axes = BuildAxes(pdfConfig_, numDimensions_);
    // put it all together and return
    BinnedED dist(name_, axes);
    const std::vector<std::string> p = pdfConfig_.GetBranchNames(numDimensions_);
    dist.SetObservables(p);
    // fill it up
    DistFiller::FillDist(dist, *data_);
    return dist;
  }

  BinnedED
  DistBuilder::BuildOscillatedDist(const std::string &name_, const int numDimensions_, const PDFConfig pdfConfig_, DataSet *data_, double deltam21_, double theta12_, std::unordered_map<int, double> indexDistance_, double& ratio)
  {
    // Create the axes
    AxisCollection axes = BuildAxes(pdfConfig_, numDimensions_);
    // Put it all together
    BinnedED dist(name_, axes);
    const std::vector<std::string> p = pdfConfig_.GetBranchNames(numDimensions_);
    dist.SetObservables(p);
    TRandom3 *rndm = new TRandom3();

    // Fill it up
    for (int i = 0; i < data_->GetNEntries(); i++)
    {
      if (!(i % 10000000))
        std::cout << i << "/" << data_->GetNEntries() << std::endl;
      Event ev = data_->GetEntry(i);
      double nuEnergy = ev.GetDatum("nu_energy");
      int reacIndex = ev.GetDatum("reactorIndex");
      double baseline = indexDistance_[reacIndex];
      double oscprob = antinufit::OscProb2(baseline, nuEnergy, deltam21_, sin(M_PI*theta12_/180)*sin(M_PI*theta12_/180));
      double r = rndm->Rndm();
      if (r > oscprob)
      {
        dist.Fill(ev);
      }
    }

    ratio = dist.Integral() / data_->GetNEntries();

    return dist;
  }
}
