#include <OscGridConfigLoader.hh>

namespace antinufit
{

  OscGridConfigLoader::OscGridConfigLoader(const std::string &filePath_)
  {
    fPath = filePath_;
  }

  OscGridConfigLoader::~OscGridConfigLoader()
  {
    ConfigLoader::Close();
  }

  OscGridConfig
  OscGridConfigLoader::Load() const
  {

    OscGridConfig ret;
    ConfigLoader::Open(fPath);
    std::string filename;
    double distance;
    double minE;
    double maxE;
    int numValsE;
    double minDm21sq;
    double maxDm21sq;
    int numValsDm21sq;
    double minSsqth12;
    double maxSsqth12;
    int numValsSsqth12;

    ConfigLoader::Load("summary", "filename", filename);
    try
    {
      ConfigLoader::Load("summary", "distance", distance);
    }
    catch (ConfigFieldMissing)
    {
      distance = 0;
    }

    ConfigLoader::Load("summary", "mine", minE);
    ConfigLoader::Load("summary", "maxe", maxE);
    ConfigLoader::Load("summary", "numvalse", numValsE);
    ConfigLoader::Load("summary", "mindm21sq", minDm21sq);
    ConfigLoader::Load("summary", "maxdm21sq", maxDm21sq);
    ConfigLoader::Load("summary", "numvalsdm21sq", numValsDm21sq);
    ConfigLoader::Load("summary", "minssqth12", minSsqth12);
    ConfigLoader::Load("summary", "maxssqth12", maxSsqth12);
    ConfigLoader::Load("summary", "numvalsssqth12", numValsSsqth12);

    ret.SetFilename(filename);
    ret.SetDistance(distance);
    ret.SetMinE(minE);
    ret.SetMaxE(maxE);
    ret.SetNumValsE(numValsE);
    ret.SetMinDm21sq(minDm21sq);
    ret.SetMaxDm21sq(maxDm21sq);
    ret.SetNumValsDm21sq(numValsDm21sq);
    ret.SetMinSsqth12(minSsqth12);
    ret.SetMaxSsqth12(maxSsqth12);
    ret.SetNumValsSsqth12(numValsSsqth12);

    return ret;
  }

}
