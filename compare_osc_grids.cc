#include <OscGrid.hh>

#include "TFile.h"
#include "TTree.h"

int main(int argc, char *argv[])
{

  if (argc != 4)
  {
    std::cout << "Usage: compare_oscgrids oscgrid_file_1 oscgrid_file_2 outfile" << std::endl;
    return 1;
  }

  std::string oscgridFile1(argv[1]);
  std::string oscgridFile2(argv[2]);
  std::string outfile(argv[3]);
  antinufit::OscGrid osc1(oscgridFile1);
  antinufit::OscGrid osc2(oscgridFile2);

  std::vector<double> fEnergyVals = osc1.GetEVector();
  std::vector<double> fDm21sqVals = osc1.GetDm21sqVector();
  std::vector<double> fSsqth12Vals = osc1.GetSsqth12Vector();
  std::vector<double> fProbVals = osc1.GetProbVector();

  TFile *f = new TFile(outfile.c_str(), "RECREATE");
  TTree *t1 = new TTree("T","Probabilities");

  double p1;
  double p2;
  double energy;
  double dm21sq;
  double ssqth21;
  double fracdiff;

  t1->Branch("p1", &p1);
  t1->Branch("p2", &p2);
  t1->Branch("energy", &energy);
  t1->Branch("dm21sq", &dm21sq);
  t1->Branch("ssqth21", &ssqth21);
  t1->Branch("fracdiff", &fracdiff);

  size_t numProbs = fProbVals.size();
  size_t i = 0;

  for (int iDmsq21 = 0; iDmsq21 < fDm21sqVals.size(); iDmsq21++)
  {
    for (int iSsqth12 = 0; iSsqth12 < fSsqth12Vals.size(); iSsqth12++)
    {
      for (int iEnergy = 0; iEnergy < fEnergyVals.size(); iEnergy++)
      {

        if (i % (numProbs / 10) == 0)
        {
          std::cout << (i * 100) / numProbs << "% events complete (" << i << "/" << numProbs << ")" << std::endl;
        }

        int probIndex = iDmsq21 * (fSsqth12Vals.size() * fEnergyVals.size()) + iSsqth12 * fEnergyVals.size() + iEnergy;

        p1 = fProbVals.at(probIndex);
        energy = fEnergyVals.at(iEnergy);
        dm21sq = fDm21sqVals.at(iDmsq21);
        ssqth21 = fSsqth12Vals.at(iSsqth12);
        p2 = osc2.Evaluate(energy, dm21sq, ssqth21);

        fracdiff = abs(p1 - p2) / p1;

        t1->Fill();

        i++;
      }
    }
  }
  f->cd();
  t1->Write();
}
