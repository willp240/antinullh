#include <BuildDist.hh>
#include <iostream>

namespace antinufit
{

  AxisCollection
  BuildDist::BuildAxes(const PDFConfig &config_)
  {

    int numAxes = config_.GetAxisCount();
    return BuildAxes(config_, numAxes);
  }

  AxisCollection
  BuildDist::BuildAxes(const PDFConfig &config_, const int numDimensions)
  {

    if (numDimensions > config_.GetAxisCount())
    {
      throw ValueError("Number of axes required is more than number declared in config");
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
  BuildDist::Build(const std::string &name_, const PDFConfig &pdfConfig_, DataSet *data_)
  {

    int numAxes = pdfConfig_.GetAxisCount();
    return Build(name_, numAxes, pdfConfig_, data_);
  }

  BinnedED
  BuildDist::Build(const std::string &name_, const int numDimensions_, const PDFConfig &pdfConfig_, DataSet *data_)
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
  BuildDist::BuildOscillatedDist(const std::string &name_, const int numDimensions_, const PDFConfig &pdfConfig_, DataSet *data_, double deltam21_, double sinsqth12_, const std::unordered_map<int, double> &indexDistance_, double &ratio)
  {
    // Create the axes
    AxisCollection axes = BuildAxes(pdfConfig_, numDimensions_);
    // Put it all together
    BinnedED dist(name_, axes);
    const std::vector<std::string> p = pdfConfig_.GetBranchNames(numDimensions_);
    dist.SetObservables(p);
    TRandom3 rndm;

    // Fill it up
    for (int i = 0; i < data_->GetNEntries(); i++)
    {
      if (!(i % 10000000))
        std::cout << i << "/" << data_->GetNEntries() << std::endl;
      Event ev = data_->GetEntry(i);
      double nuEnergy = ev.GetDatum("nu_energy");
      int reacIndex = ev.GetDatum("reactorIndex");
      double baseline = indexDistance_.at(reacIndex);
      double oscprob = antinufit::OscProb(baseline, nuEnergy, deltam21_, sinsqth12_);
      double r = rndm.Rndm();
      if (r > oscprob)
      {
        dist.Fill(ev);
      }
    }

    ratio = dist.Integral() / data_->GetNEntries();

    return dist;
  }

  BinnedED
  BuildDist::BuildFlatDist(const std::string &name_, const int numDimensions_, const PDFConfig &pdfConfig_)
  {
    // Create the axes
    AxisCollection axes = BuildAxes(pdfConfig_, numDimensions_);
    // put it all together and return
    BinnedED dist(name_, axes);
    const std::vector<std::string> p = pdfConfig_.GetBranchNames(numDimensions_);
    dist.SetObservables(p);

    // fill it up
    for (size_t ibin = 0; ibin < dist.GetNBins(); ibin++)
    {
      dist.SetBinContent(ibin, 1 / dist.GetNBins());
    }

    return dist;
  }
}
