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

void Po214alphagammaFit(){
    int simbins = 38; double startbin = 0.6; double endbin =2.5;

    TCanvas *c1 = new TCanvas("c1","",800,500);
    //c1->SetLogy();

    std::vector<double> mcR; std::vector<double> R ; std::vector<double> R_4m; std::vector<std::string> runid;
    std::string runlist; std::string Eventtype; 
    //runlist =  "/home/huangp/AntiNu/preliminary_goldrunlist.txt";
    //runlist =  "/home/huangp/AntiNu/antinu_runlist_UPDATED_copy.txt";
    runlist =  "/home/huangp/RafDownload/DownloadList/2p2Po214_list.dat";
    
    runid = readrunidfromtable(runlist);

    TChain* Tc = new TChain("output"); 
    Int_t evIndex;
    double  energy,delayedecorr ;
    Bool_t fitValid;
    Tc->SetBranchAddress("fitValid", &fitValid);
    Tc->SetBranchAddress("evIndex", &evIndex);
    Tc->SetBranchAddress("energy", &energy);
    TNtuple* nt = new TNtuple("nt","nt","delayedEcorr"); 
    std::string filepath =  "/data/snoplus2/weiiiii/2p2Po214/Ntuple/AlphaGammaSim/";
    for (size_t irun = 0; irun < runid.size(); ++irun){
        std::string input_file  = "Po214*"+runid[irun]+"*.root";

        //std::cout<<(filepath+input_file).c_str()<<std::endl;
        
        Tc->Add( (filepath+input_file).c_str());
        
        // check this run number exists or not
        if(Tc->GetEntries() == 0){
            std::cerr << "No such file or Error opening recon file in run: "<< runid[irun]<< std::endl;
            //exit(-1);
        }
    }
    //std::cout<<"Entries "<<Tc->GetEntries()<<std::endl;
    TH1D *Data_energyh1  = new TH1D("Data_energyh1","",simbins,startbin ,endbin);
    for (int iEntry = 0; iEntry < Tc->GetEntries(); iEntry++){
        Tc->GetEntry(iEntry);
        //Data_energyh1->Fill(delayedEcorr);
        if(evIndex==0 && fitValid!=0 && energy>1.2 ){
            delayedecorr = energy;
            nt->Fill(delayedecorr);
        }
    }
    nt->Draw("delayedEcorr");
    // ROOFIT
    
    RooRealVar delayedEcorr("delayedEcorr","delayedEcorr",startbin,endbin);

    //RooDataHist energybindata("energybindata","energybindata",RooArgList(delayedEcorr),Data_energyh1);
    
    RooDataSet data("data","data",nt,RooArgSet(delayedEcorr)); 
    data.Print("v"); 

    // Construct Gaus PDF

    RooRealVar mu_alphagamma("mu_alphagamma", "mean parameter", 1.3, 1.45, 1.6);
    RooRealVar sigma_alphagamma("sigma_alphagamma", "width parameter", 0.5, 0.05, 1.0);

    RooGaussian gaus_alphagamma("gaus_alphagamma", "214 alphagamma Gaussian", delayedEcorr, mu_alphagamma, sigma_alphagamma); 
    gaus_alphagamma.fitTo(data);
    RooPlot* frame = delayedEcorr.frame();
    data.plotOn(frame,Name("data"),Binning(38)) ;
    gaus_alphagamma.plotOn(frame,Name("gaus_alphagamma"),Binning(38)) ;
    frame->Draw("E0");
    TLegend *leg1 = new TLegend(0.5, 0.68, 0.90, 0.9);

    leg1->SetTextSize(.04);
    leg1->SetTextFont(42);
    //leg1->SetFillColor(kWhite);
    //leg1->SetLineColor(kWhite);
    
    leg1->AddEntry(frame->findObject("data"),"Data: Po214 alpha+gamma","P");
    leg1->AddEntry(frame->findObject("gaus_alphagamma"),"214 alphagamma Gaussian", "L");
    leg1->Draw("same");  
    
       //save into workspace
    RooWorkspace *w_alphagamma = new RooWorkspace("w_alphagamma", "workspace");
 
    // Import model and all its components into the workspace
    w_alphagamma->import(gaus_alphagamma);
 
    // Import data into the workspace
    w_alphagamma->import(data);
 
    // Print workspace contents
    w_alphagamma->Print();

    // S a v e   w o r k s p a c e   i n   f i l e
   // -------------------------------------------
 
   // Save the workspace into a ROOT file
   w_alphagamma->writeToFile("Po214alphagammaFit.root");
 
   delete w_alphagamma;
    

}


