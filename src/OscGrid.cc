#include "OscGrid.hh"

namespace antinufit
{

  OscGrid::OscGrid(const std::string &filename_, const double &distance_)
  {

    fFilename = filename_;
    fDistance = distance_;
    Load();
  }

  OscGrid::Load()
  {

    std::cout << "Loading probability information from " << fFilename << std::endl;
    std::ifstream infile;
    infile.open(fFilename, std::ios_base::in);

    // For each axis, load the scan points
    fEnergyVals.clear();
    fTheta12Vals.clear();
    fDmsq21Vals.clear();

    std::vector<std::vector<double> *> axes = {&fEnergyVals, &fTheta12Vals, &fDmsq21Vals};

    for (auto &axis : axes)
    {
      std::string line;
      std::getline(infile, line);
      std::istringstream line_s(line);
      while (line_s)
      {
        double val;
        line_s >> val;
        axis->push_back(val);
      }
      axis->pop_back();
    }
    fMinE = fEnergyVals.at(0);
    fMaxE = fEnergyVals.back();
    fNumValsE = fEnergyVals.size();
    fMinDm21sq = fDmsq21Vals.at(0);
    fMaxDm21sq = fDmsq21Vals.back();
    fNumValsDm21sq = fDmsq21Vals.size();
    fMinSsqth12 = fTheta12Vals.at(0);
    fMaxSsqth12 = fTheta12Vals.back();
    fNumValsSsqth12 = fTheta12Vals.size();

    // Load the probability values
    for (std::string line; std::getline(infile, line);)
    {
      fProbVals.push_back(std::stod(line));
    }
    std::cout << "Loading of OscGrid complete!" << std::endl << std::endl;
  }

  void OscGrid::CalcGrid()
  {

    fEnergyVals = linspace(fMinE, fMaxE, fNumValsE);
    fDm21sqVals = linspace(fMinDm21sq, fMaxDm21sq, fNumValsDm21sq);
    fSsqth12Vals = linspace(fMinSsqth12, fMaxSsqth12, fNumValsSsqth12);

    const size_t numVals = fNumValsE * fNumValsDm21sq * fNumValsSsqth12;
    fProbVals.reserve(numVals);
    size_t i;

    for (const auto &dmsq21 : fDm21sqVals)
    {
      for (const auto &theta12 : fSsqth12Vals)
      {
        for (const auto &energy : fEnergyVals)
        {
          if (i % (numVals / 100) == 0)
          {
            std::cout << (i * 100) / numVals << "% events complete (" << i << "/" << numVals << ")\n";
          }

          double prob = OscProb(distance, energy, dmsq21, theta12)
                            fProb_vals.push_back(prob);
          i++;
        }
      }
    }
  }

  void OscGrid::Write()
  {

    std::cout << "Writing Oscillation Grid information to " << fFilename << std::endl;
    // Open output file for writing
    std::ofstream outfile;
    outfile.open(fFilename, std::ios_base::out);
    // For each axis, write the scan points
    const auto axes = {fEnergyVals, fDm21sqVals, fSsqth12Vals};
    for (const auto &axis : axes)
    {
      for (const auto &val : axis)
      {
        outfile << val << " ";
      }
      outfile << std::endl;
    }

    // Now write the probability vals
    for (const auto &prob : fProbVals)
    {
      outfile << prob << std::endl;
    }
    std::cout << "Writing complete!" << std::end;
  }

  double OscGrid::Evaluate()
  {
  }
}
