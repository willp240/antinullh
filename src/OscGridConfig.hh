#ifndef __ANTINUFIT__OscGridConfig__
#define __ANTINUFIT__OscGridConfig__

// OXO headers
#include <ContainerTools.hpp>

namespace antinufit
{
  class OscGridConfig
  {
  public:
    void SetFilename(std::string);
    std::string GetFilename() const;
    void SetDistance(double);
    double GetDistance() const;
    void SetMinE(double);
    double GetMinE() const;
    void SetMaxE(double);
    double GetMaxE() const;
    void SetNumValsE(int);
    int GetNumValsE() const;
    void SetMinDm21sq(double);
    double GetMinDm21sq() const;
    void SetMaxDm21sq(double);
    double GetMaxDm21sq() const;
    void SetNumValsDm21sq(int);
    int GetNumValsDm21sq() const;
    void SetMinSsqth12(double);
    double GetMinSsqth12() const;
    void SetMaxSsqth12(double);
    double GetMaxSsqth12() const;
    void SetNumValsSsqth12(int);
    int GetNumValsSsqth12() const;

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
  };
}
#endif
