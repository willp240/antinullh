#include <FitSummary.hh>
#include <algorithm>
#include <iomanip>
#include <iostream>

namespace antinufit
{
  void PrintFitSummary(const FitParameters &p,
                       const std::map<std::string, std::vector<std::string>> &datasets,
                       const OscParams &osc,
                       const PDFConfig &pdfConfig)
  {

    // Preferred ordering for readability
    std::vector<std::string> orderedNames = {
        "deltam21",
        osc.theta12name,
        "reactor_nubar",
        "geonu_U",
        "geonu_Th",
        "alphan_CScatter",
        "alphan_OExcited",
        "alphan_PRecoil",
        "sideband",
        "energy_scale",
        "energy_conv",
        "birks_constant",
        "p_recoil_energy_scale"};

    // Append any remaining params not in the preferred list
    for (ParameterDict::const_iterator it = p.noms.begin(); it != p.noms.end(); ++it)
    {
      if (std::find(orderedNames.begin(), orderedNames.end(), it->first) == orderedNames.end())
        orderedNames.push_back(it->first);
    }

    std::cout << "\n";
    std::cout << "************** Fit Parameters **************\n";
    std::cout << " -------------------------------------------------------------------------------------------------------------------------------------------------------\n";
    std::cout << "| " << std::left << std::setw(25) << "Name"
              << "| " << std::left << std::setw(15) << "Nominal"
              << "| " << std::left << std::setw(15) << "Minimum"
              << "| " << std::left << std::setw(15) << "Maximum"
              << "| " << std::left << std::setw(20) << "Datasets"
              << "| " << std::left << std::setw(20) << "Constraint Mean"
              << "| " << std::left << std::setw(20) << "Constraint Sigma"
              << "| " << std::left << std::setw(6) << "Fixed"
              << "| \n";
    std::cout << " =======================================================================================================================================================\n";

    for (std::vector<std::string>::const_reference name : orderedNames)
    {
      if (p.noms.find(name) == p.noms.end())
        continue;

      std::cout << "| " << std::left << std::setw(25) << name
                << "| " << std::left << std::setw(15) << p.noms.at(name)
                << "| " << std::left << std::setw(15) << p.mins.at(name)
                << "| " << std::left << std::setw(15) << p.maxs.at(name)
                << "| ";

      // Datasets column
      std::string datasetString;
      const std::map<std::string, std::vector<std::string>>::const_iterator dsIt = datasets.find(name);
      if (dsIt != datasets.end())
        datasetString = Join(dsIt->second);
      std::cout << std::left << std::setw(20) << datasetString << "| ";

      // Constraint columns (simple Gaussian priors only)
      const ParameterDict::const_iterator cmIt = p.constraintMeans.find(name);
      if (cmIt != p.constraintMeans.end())
      {
        std::cout << std::left << std::setw(20) << cmIt->second
                  << "| " << std::left << std::setw(20) << p.constraintSigmas.at(name)
                  << "| ";
      }
      else
      {
        std::cout << std::left << std::setw(20) << ""
                  << "| " << std::left << std::setw(20) << ""
                  << "| ";
      }

      bool fixed = false;
      const std::map<std::string, bool>::const_iterator fxIt = p.fixedPars.find(name);
      if (fxIt != p.fixedPars.end())
        fixed = fxIt->second;

      std::cout << std::left << std::setw(6) << fixed << "| \n";
    }

    std::cout << " -------------------------------------------------------------------------------------------------------------------------------------------------------\n";

    // Ratio constraints
    if (!p.ratioMeans.empty())
    {
      std::cout << "\n************** Ratio Constraints **************\n";
      std::cout << " -------------------------------------------------------------------------------------------------\n";
      std::cout << "| " << std::left << std::setw(25) << "Parameter"
                << "| " << std::left << std::setw(25) << "Parameter"
                << "| " << std::left << std::setw(20) << "Constraint Mean"
                << "| " << std::left << std::setw(20) << "Constraint Sigma"
                << "| \n";
      std::cout << " =================================================================================================\n";

      for (ParameterDict::const_iterator it = p.ratioMeans.begin(); it != p.ratioMeans.end(); ++it)
      {
        const std::string &a = it->first;
        std::cout << "| " << std::left << std::setw(25) << a
                  << "| " << std::left << std::setw(25) << p.ratioParName.at(a)
                  << "| " << std::left << std::setw(20) << it->second
                  << "| " << std::left << std::setw(20) << p.ratioSigmas.at(a)
                  << "| \n";
      }
      std::cout << " -------------------------------------------------------------------------------------------------\n";
    }

    // Fraction constraints
    if (!p.fracMeans.empty())
    {
      std::cout << "\n************** Fractional Diff. Constraints **************\n";
      std::cout << " -------------------------------------------------------------------------------------------------\n";
      std::cout << "| " << std::left << std::setw(25) << "Parameter"
                << "| " << std::left << std::setw(25) << "Parameter"
                << "| " << std::left << std::setw(20) << "Constraint Mean"
                << "| " << std::left << std::setw(20) << "Constraint Sigma"
                << "| \n";
      std::cout << " =================================================================================================\n";

      for (ParameterDict::const_iterator it = p.fracMeans.begin(); it != p.fracMeans.end(); ++it)
      {
        const std::string &a = it->first;
        std::cout << "| " << std::left << std::setw(25) << a
                  << "| " << std::left << std::setw(25) << p.fracParName.at(a)
                  << "| " << std::left << std::setw(20) << it->second
                  << "| " << std::left << std::setw(20) << p.fracSigmas.at(a)
                  << "| \n";
      }
      std::cout << " -------------------------------------------------------------------------------------------------\n";
    }

    // Shape constraints
    if (!p.shapeMeans.empty())
    {
      std::cout << "\n************** Shape Constraints **************\n";
      std::cout << " -------------------------------------------------------------------------------------------------\n";
      std::cout << "| " << std::left << std::setw(25) << "Key"
                << "| " << std::left << std::setw(25) << "Shape Function"
                << "| " << std::left << std::setw(20) << "Constraint Mean"
                << "| " << std::left << std::setw(20) << "Constraint Sigma"
                << "| \n";
      std::cout << " =================================================================================================\n";

      for (ParameterDict::const_iterator it = p.shapeMeans.begin(); it != p.shapeMeans.end(); ++it)
      {
        const std::string &k = it->first;
        std::cout << "| " << std::left << std::setw(25) << k
                  << "| " << std::left << std::setw(25) << p.shapeFuncName.at(k)
                  << "| " << std::left << std::setw(20) << it->second
                  << "| " << std::left << std::setw(20) << p.shapeSigmas.at(k)
                  << "| \n";

        const std::map<std::string, std::vector<std::string>>::const_iterator pnIt = p.shapeParNames.find(k);
        if (pnIt != p.shapeParNames.end())
        {
          for (std::vector<std::string>::const_reference pname : pnIt->second)
          {
            std::cout << "| " << std::left << std::setw(25) << ("  " + pname)
                      << "| " << std::left << std::setw(25) << ""
                      << "| " << std::left << std::setw(20) << ""
                      << "| " << std::left << std::setw(20) << ""
                      << "| \n";
          }
        }
      }
      std::cout << " -------------------------------------------------------------------------------------------------\n";
    }

    // Correlations
    if (!p.corrs.empty())
    {
      std::cout << "\n************** Correlations **************\n";
      std::cout << " ---------------------------------------------------------------------------\n";
      std::cout << "| " << std::left << std::setw(25) << "Parameter"
                << "| " << std::left << std::setw(25) << "Parameter"
                << "| " << std::left << std::setw(20) << "Correlation"
                << "| \n";
      std::cout << " ===========================================================================\n";

      for (ParameterDict::const_iterator it = p.corrs.begin(); it != p.corrs.end(); ++it)
      {
        const std::string &a = it->first;
        std::cout << "| " << std::left << std::setw(25) << a
                  << "| " << std::left << std::setw(25) << p.corrParName.at(a)
                  << "| " << std::left << std::setw(20) << it->second
                  << "| \n";
      }
      std::cout << " ---------------------------------------------------------------------------\n";
    }

    std::cout << "\nLLH buffer bins by axis:\n";
    for (const std::string &axisName : pdfConfig.GetAxisNames())
    {
      if (pdfConfig.HasLLHBufferBins(axisName))
      {
        const std::pair<int, int> bins = pdfConfig.GetLLHBufferBins(axisName);
        std::cout << "  " << axisName << ": low=" << bins.first
                  << ", high=" << bins.second << '\n';
      }
      else
        std::cout << "  " << axisName << ": none\n";
    }
    std::cout << "\n";
  }
} // namespace antinufit
