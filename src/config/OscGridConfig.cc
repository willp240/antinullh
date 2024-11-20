#include "OscGridConfig.hh"

namespace antinufit
{

  std::string OscGridConfig::GetFilename() const
  {
    return fFilename;
  }

  void
  OscGridConfig::SetFilename(std::string filename_)
  {
    fFilename = filename_;
  }

    std::string OscGridConfig::GetReactorsJsonFile() const
  {
    return fReactorsJsonFilename;
  }

  void
  OscGridConfig::SetReactorsJsonFile(std::string reactorsjsonfilename_)
  {
    fReactorsJsonFilename = reactorsjsonfilename_;
  }

  double OscGridConfig::GetMinE() const
  {
    return fMinE;
  }

  void
  OscGridConfig::SetMinE(double mine_)
  {
    fMinE = mine_;
  }

  double OscGridConfig::GetMaxE() const
  {
    return fMaxE;
  }

  void
  OscGridConfig::SetMaxE(double maxe_)
  {
    fMaxE = maxe_;
  }

  int OscGridConfig::GetNumValsE() const
  {
    return fNumValsE;
  }

  void
  OscGridConfig::SetNumValsE(int numvalse_)
  {
    fNumValsE = numvalse_;
  }

  double OscGridConfig::GetMinDm21sq() const
  {
    return fMinDm21sq;
  }

  void
  OscGridConfig::SetMinDm21sq(double minDm21sq_)
  {
    fMinDm21sq = minDm21sq_;
  }

  double OscGridConfig::GetMaxDm21sq() const
  {
    return fMaxDm21sq;
  }

  void
  OscGridConfig::SetMaxDm21sq(double maxdm21sq_)
  {
    fMaxDm21sq = maxdm21sq_;
  }

  int OscGridConfig::GetNumValsDm21sq() const
  {
    return fNumValsDm21sq;
  }

  void
  OscGridConfig::SetNumValsDm21sq(int numvalsdm21sq_)
  {
    fNumValsDm21sq = numvalsdm21sq_;
  }

  double OscGridConfig::GetMinSsqth12() const
  {
    return fMinSsqth12;
  }

  void
  OscGridConfig::SetMinSsqth12(double minSsqth12_)
  {
    fMinSsqth12 = minSsqth12_;
  }

  double OscGridConfig::GetMaxSsqth12() const
  {
    return fMaxSsqth12;
  }

  void
  OscGridConfig::SetMaxSsqth12(double maxSsqth12_)
  {
    fMaxSsqth12 = maxSsqth12_;
  }

  int OscGridConfig::GetNumValsSsqth12() const
  {
    return fNumValsSsqth12;
  }

  void
  OscGridConfig::SetNumValsSsqth12(int numvalsSsqth12_)
  {
    fNumValsSsqth12 = numvalsSsqth12_;
  }

}
