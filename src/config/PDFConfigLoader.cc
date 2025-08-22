#include <PDFConfigLoader.hh>

namespace antinufit
{

  PDFConfigLoader::PDFConfigLoader(const std::string &filePath_)
  {
    fPath = filePath_;
  }

  PDFConfig
  PDFConfigLoader::Load() const
  {
    ConfigLoader::Open(fPath);

    std::set<std::string> toLoad = ConfigLoader::ListSections();
    toLoad.erase("summary");

    std::vector<std::string> order;
    ConfigLoader::Load("summary", "build_order", order);

    std::vector<std::string> data_axes;
    ConfigLoader::Load("summary", "data_axes", data_axes);

    PDFConfig retVal;

    double min;
    double max;
    std::string name;
    std::string branchName;
    std::string texName;
    int binCount;
    for (size_t i = 0; i < order.size(); i++)
    {
      if (std::find(toLoad.begin(), toLoad.end(), order.at(i)) == toLoad.end())
        throw NotFoundError(Formatter() << "PDFConfigLoader:: " << order.at(i)
                                        << " is in the build order but has no section!");

      name = order.at(i);
      ConfigLoader::Load(name, "min", min);
      ConfigLoader::Load(name, "max", max);
      ConfigLoader::Load(name, "n_bins", binCount);
      ConfigLoader::Load(name, "branch_name", branchName);
      ConfigLoader::Load(name, "tex_name", texName);

      // Remove "" if it surrounds the tex name string from the config
      texName = stripQuoteMarks(texName);

      retVal.AddAxis(name, branchName, texName, binCount, min, max);
    }
    retVal.SetDataBranchNames(data_axes);
    return retVal;
  }

  PDFConfigLoader::~PDFConfigLoader()
  {
    ConfigLoader::Close();
  }

}
