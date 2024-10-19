#ifndef __ANTINUFIT__OscGrid__
#define __ANTINUFIT__OscGrid__

#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>

#include "linearinterpolation.hxx"

#include <Utilities.hh>
#include <Functions.hh>

#include <TH3D.h>

namespace antinufit
{
  class OscGrid
  {
  public:
    OscGrid(const std::string &filename_, const double &distance_, const double &mine_, const double &maxe_, const double &numvalse_,
            const double &mindm21sq_, const double &maxdm21sq_, const double &numvalsdm21sq_,
            const double &minssqth12_, const double &maxssqth12_, const double &numvalsssqth12_)
        : fFilename(filename_), fDistance(distance_), fMinE(mine_), fMaxE(maxe_), fNumValsE(numvalse_),
          fMinDm21sq(mindm21sq_), fMaxDm21sq(maxdm21sq_), fNumValsDm21sq(numvalsdm21sq_),
          fMinSsqth12(minssqth12_), fMaxSsqth12(maxssqth12_), fNumValsSsqth12(numvalsssqth12_) {};

    OscGrid(const std::string &filename_, const double &distance_);
    OscGrid() {};

    void Load();
    void CalcGrid();
    void Write();
    double Evaluate(double, double, double);
    TH3D *MakeHist();

  private:
    std::string fFilename;
    double fDistance;

    double fMinE;
    double fMaxE;
    int fNumValsE;
    double fMinDm21sq;
    double fMaxDm21sq;
    int fNumValsDm21sq;
    double fMinSsqth12;
    double fMaxSsqth12;
    int fNumValsSsqth12;

    std::vector<double> fEnergyVals;
    std::vector<double> fDm21sqVals;
    std::vector<double> fSsqth12Vals;
    std::vector<double> fProbVals;

    TH3D *hist;
  };
}
#endif
