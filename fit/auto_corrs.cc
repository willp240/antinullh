#include <iostream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <TMath.h>
#include <string>

#include "TObjArray.h"
#include "TChain.h"
#include "TFile.h"
#include "TBranch.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TLegend.h"
#include "TString.h"
#include "TStyle.h"
#include "TH1.h"
#include "TF1.h"
#include "TH2.h"
#include "TGraphErrors.h"
#include "TVectorD.h"
#include "TColor.h"
#include "TKey.h"
#include "TPaveText.h"

#include <FitConfigLoader.hh>
#include <FitConfig.hh>

using namespace bbfit;

int main(int argc, char *argv[]) {
  
  if (argc != 2) {
    std::cerr << "./bin/auto_corrs root_file_to_analyse.root " << std::endl;
    exit(-1);
  }

  std::string fileName = argv[1];
  TFile* inputFile = new TFile(fileName.c_str(), "OPEN");

  std::string outFileName = fileName.substr(0, fileName.find(".root"));;
  outFileName+= "_autoCorr.root";

  TChain* pos = new TChain("posteriors", "");
  pos->Add(fileName.c_str());

  TObjArray* branchList = (TObjArray*)pos->GetListOfBranches();
  int nBranches = branchList->GetEntries();
  int nSteps = -999;
  int nPars = 0;
  int maxLag = 1000;
  std::vector<std::string> parNames;

  for (int i = 0; i < nBranches; i++) {
    TBranch *branch = (TBranch*)branchList->At(i);
    std::string branchName = std::string(branch->GetName());
    
    if (nSteps == -999)
      nSteps = branch->GetEntries();
    
    if(branchName!= "Step" && branchName != "StepTime" && branchName != "Accepted" && branchName != "LogL"){
      nPars++;
      parNames.push_back(branchName);
    }
  }

  double **parVals = new double*[nSteps]();
  double *parSums = new double[nPars]();;

  for (int i = 0; i < nSteps; ++i) {
    parVals[i] = new double[nPars]();
    for (int j = 0; j < nPars; ++j) {
      parVals[i][j] = -999.99;
      parSums[j] = 0.0;
    }
  }  

  pos->SetBranchStatus("*", false);
  for (int i = 0; i < nPars; ++i) {
    pos->SetBranchStatus(parNames[i].c_str(), true);
  }

  // Loop over all steps
  for (int i = 0; i < nSteps; ++i) {
    if (i % 100000 == 0) {
      std::cout << "Read " << i << " of " << nSteps << std::endl;
    }

    for (int j = 0; j < nPars; ++j) {
      pos->SetBranchAddress(parNames[j].c_str(), &parVals[i][j]);
    }

    pos->GetEntry(i);
    for (int j = 0; j < nPars; ++j)
      parSums[j] += parVals[i][j];
  }

  for (int i = 0; i < nPars; ++i)
    parSums[i] /= nSteps;

  delete pos;

  TFile* outputFile = new TFile(outFileName.c_str(), "RECREATE");

  double **denomSum = new double*[nPars]();
  double **numSum = new double*[nPars]();
  for (int i = 0; i < nPars; ++i) {
    denomSum[i] = new double[maxLag];
    numSum[i] = new double[maxLag];
  }

  TH1D** Lags = new TH1D*[nPars];

  for (int j = 0; j < nPars; j++) {
    for (int k = 0; k < maxLag; k++) {
      numSum[j][k] = 0.0;
      denomSum[j][k] = 0.0;
    }

    Lags[j] = new TH1D(parNames[j].c_str(), parNames[j].c_str(), maxLag, 0.0, maxLag);
    Lags[j]->GetXaxis()->SetTitle("Lag");
    Lags[j]->GetYaxis()->SetTitle("Auto-correlation");
  }

  for (int k = 0; k < maxLag; ++k) {
    for (int i = 0; i < nSteps; ++i) {
      for (int j = 0; j < nPars; ++j) {

        double Diff = parVals[i][j]-parSums[j];

        if (i < nSteps-k)
          numSum[j][k] += Diff*(parVals[i+k][j]-parSums[j]);

        denomSum[j][k] += Diff*Diff;
      }
    }
  }

  outputFile->cd();
  TDirectory *AutoCorrs = outputFile->mkdir("Auto_Corrs");
  for (int j = 0; j < nPars; ++j) {
    for (int k = 0; k < maxLag; ++k)
      Lags[j]->SetBinContent(k, numSum[j][k]/denomSum[j][k]);
    AutoCorrs->cd();
    Lags[j]->Write();
    delete Lags[j];
  }
  delete[] Lags;

  for (int i = 0; i < nSteps; ++i) {
    delete parVals[i];
  }
  delete[] parVals;

  delete parSums;

  return 0;
}

