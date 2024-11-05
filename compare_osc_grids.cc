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

  double p1;
  double p2;
  double energy;
  double dm21sq;
  double ssqth21;
  double fracdiff;

  std::vector<double> fEnergyVals = osc1.GetEVector();
  std::vector<double> fDm21sqVals = osc1.GetDm21sqVector();
  std::vector<double> fSsqth12Vals = osc1.GetSsqth12Vector();
  TTree *fProbTree = osc1.GetProbTree();
  fProbTree->SetBranchAddress("prob", &p1);
  fProbTree->SetBranchAddress("energy", &energy);
  fProbTree->SetBranchAddress("dmsq21", &dm21sq);
  fProbTree->SetBranchAddress("ssqth12", &ssqth21);

  TFile *f = new TFile(outfile.c_str(), "RECREATE");
  TTree *t1 = new TTree("T", "Probabilities");

  t1->Branch("p1", &p1);
  t1->Branch("p2", &p2);
  t1->Branch("energy", &energy);
  t1->Branch("dm21sq", &dm21sq);
  t1->Branch("ssqth21", &ssqth21);
  t1->Branch("fracdiff", &fracdiff);

  size_t numProbs = fProbTree->GetEntries();
  size_t i = 0;

  for (int iEntry = 0; iEntry < numProbs; iEntry++)
  {

    //if (iEntry % (numProbs / 1) == 0)
    {
      double frac = 100 * (double)iEntry / numProbs;
      std::cout << frac << "% events complete (" << iEntry << "/" << numProbs << ")" << std::endl;
    }

    fProbTree->GetEntry(iEntry);
    p2 = osc2.Evaluate(energy, dm21sq, ssqth21);
    fracdiff = abs(p1 - p2) / p1;
    f->cd();
    t1->Fill();
  }

  f->cd();
  t1->Write();
}
