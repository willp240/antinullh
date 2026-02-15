// usage: root -l 'effcalculation.cpp(7015,214)' or root -l 'effcalculation.cpp(7015,212)'
// Do the efficiency calculation of reactor ibd, geoantinu and (alpha,n)
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TRandom.h"
#include "TTree.h"
#include "TChain.h"
#include "TMath.h"

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <glob.h>
#include <RAT/DU/DSReader.hh>                                               
#include <RAT/DU/Utility.hh>
#include <RAT/DS/Entry.hh>
#include <RAT/GeoUtils.hh>

#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooConstVar.h"
#include "RooPolynomial.h"
#include "RooKeysPdf.h"
#include "RooNDKeysPdf.h"
#include "RooProdPdf.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TH1.h"
#include "RooPlot.h"
#include "RooWorkspace.h"
using namespace RooFit ;

std::vector<std::string>  readrunidfromtable(const std::string& filename){
    std::vector<std::string> indexes;
    //std::ifstream inputFile("/home/huangp/BiPo/PhysicsRunList/bismsb_bipo214_simlist.txt");
    std::ifstream inputFile(filename);

    if (!inputFile.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return indexes;
    }

    std::string index;
    while (inputFile >> index) {
        //std::cout << index <<std::endl;
        indexes.push_back(index);
    }

    inputFile.close();
    return indexes;

} 

void sideband_copy(){
    TCanvas *c1 = new TCanvas("c1","",800,500);
    c1->SetLogy();

    std::vector<double> mcR; std::vector<double> R ; std::vector<double> R_4m; std::vector<std::string> runid; std::vector<std::string> partialrunid;
    std::string runlist; std::string Eventtype;  std::string partialrunlist;
    runlist        =       "/home/huangp/AntiNu/preliminary_goldrunlist.txt";
    partialrunlist =       "/home/huangp/AntiNu/preliminary_partialantinu.txt";
    
    runid = readrunidfromtable(runlist);
    partialrunid = readrunidfromtable(partialrunlist);
    TChain* Tc = new TChain("PoT"); 
    double Po_R, delayedecorr,Po_Energy;
    Tc->SetBranchAddress("Po_R", &Po_R);
    Tc->SetBranchAddress("Po_energy", &Po_Energy);
    TNtuple* nt = new TNtuple("nt","nt","delayedEcorr"); 
    std::string filepath        =  "/data/snoplus3/weiii/RnPo215/fullFill/rat-7.0.8/Ntuple/";
    std::string partialfilepath =  "/data/snoplus3/weiii/RnPo215/partialFill/rat-6.18.9/Ntuple/";
    
    for (size_t irun = 0; irun < runid.size(); ++irun){
        std::string input_file  = runid[irun]+"*.ntuple.root";
        Tc->Add( (filepath+input_file).c_str());
    }
    
    for (size_t irun = 0; irun < partialrunid.size(); ++irun){
        std::string input_file  = partialrunid[irun]+"*.ntuple.root";
        Tc->Add( (partialfilepath+input_file).c_str());
    }
    if(Tc->GetEntries() == 0){
            std::cerr << "Error Opening File "<<std::endl;
            //exit(-1);
        
    }
    
    //TH1D *Data_energyh1  = new TH1D("Data_energyh1","",100,0.6,2.5);
    for (int iEntry = 0; iEntry < Tc->GetEntries(); iEntry++){
        Tc->GetEntry(iEntry);
        //Data_energyh1->Fill(Po_energy);
        if(Po_R<5000){
            delayedecorr = Po_Energy;
            nt->Fill(delayedecorr);
        }
    }
    nt->Draw("delayedEcorr");
    /*
    TTree* newTree = tree->CloneTree(0);  // Clone structure only to new branch delayedEcorr

    // Define a new branch
    double delayedEcorr;
    newTree->Branch("delayedEcorr", &delayedEcorr);

    // Set branch address for the original tree
    double energy;
    nt->SetBranchAddress("Energy", &energy);

    // Loop over entries and fill the new tree
    Long64_t nentries = tree->GetEntries();
    for (Long64_t i = 0; i < nentries; i++) {
        tree->GetEntry(i);
        delayedEcorr = energy;  // Copy value
        newTree->Fill();
    }
    // ROOFIT
    */

    RooRealVar delayedEcorr("delayedEcorr","delayedEcorr",0.6,2.5);
    //RooRealVar Po_energy("Po_energy","Po_energy",0.6,2.5);
    delayedEcorr.setBins(50);

    //RooDataHist energybindata("energybindata","energybindata",RooArgList(Po_energy),Data_energyh1);
    RooDataSet data("data","data",nt,RooArgSet(delayedEcorr)); 
    data.Print("v"); 

    // Construct KDE PDF
    RooKeysPdf k1("k1","k1",delayedEcorr,data,RooKeysPdf::NoMirror) ;

    RooPlot* frame = delayedEcorr.frame();
    //frame->SetMinimum(0.1);
    data.plotOn(frame,Name("data")) ;
    k1.plotOn(frame,Name("k1")) ;
    frame->Draw("E0");
    TLegend *leg1 = new TLegend(0.5, 0.68, 0.90, 0.9);

    leg1->SetTextSize(.04);
    leg1->SetTextFont(42);
    //leg1->SetFillColor(kWhite);
    //leg1->SetLineColor(kWhite);
    
    leg1->AddEntry(frame->findObject("data"),"Data","P");
    leg1->AddEntry(frame->findObject("k1"),"Kernal Density Estimation", "L");
    leg1->Draw("same");  

    
    
    

    
    //save into workspace
    RooWorkspace *w_Po215 = new RooWorkspace("w_Po215", "workspace");
 
    // Import model and all its components into the workspace
    w_Po215->import(k1);
 
    // Import data into the workspace
    w_Po215->import(data);
 
    // Print workspace contents
    w_Po215->Print();

    // S a v e   w o r k s p a c e   i n   f i l e
   // -------------------------------------------
 
   // Save the workspace into a ROOT file
   w_Po215->writeToFile("sidebandkdefit.root");
 
   delete w_Po215;
    
   double x_min = 1.735;
   double x_max = 2.5;

   // Set the range for x variable
   delayedEcorr.setRange("region", x_min, x_max);

   // Now, to integrate a particular component (for example, gaus1), do the following:
   
   // First, create an integral of the specific component (gaus1) over the region
   RooAbsReal* integral_k1 = k1.createIntegral(RooArgSet(delayedEcorr), Range("region"));
   
   // Get the integral value (total probability in the region for gaus1)
   double integratedValue_k1 = integral_k1->getVal();

   // Print the integrated value of gaus1
   std::cout << "Integrated value of KDE in the region [" << x_min << ", " << x_max << "] is: " 
            << integratedValue_k1 << std::endl;
    std::cout << "Integrated value of KDE outside the region [" << x_min << ", " << x_max << "] is: " 
            << 1 - integratedValue_k1 << std::endl;
}


