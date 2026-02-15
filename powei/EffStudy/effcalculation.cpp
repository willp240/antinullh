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

std::vector<std::string> glob(const char *pattern) {
    glob_t g;
    glob(pattern, GLOB_TILDE, nullptr, &g); // one should ensure glob returns 0!
    std::vector<std::string> filelist;
    filelist.reserve(g.gl_pathc);// reserve enough space for putting each filename inside vector
    for (size_t i = 0; i < g.gl_pathc; ++i) {
        filelist.emplace_back(g.gl_pathv[i]);//append file in vector
    }
    return filelist;
}
double GetZAVoffset(RAT::DB *db,RAT::DS::Run &run,int runID){
    
    
    run.SetRunID(runID);
    db->BeginOfRun(run);
    std::vector<double> AVPos = RAT::GeoUtil::UpdateAVOffsetVectorFromDB();
    double zOff = AVPos[2];
    std::cout << "AV Shift is: " << zOff << std::endl;
    return zOff;
}
void GetZFromOriginal(RAT::DB *db,RAT::DS::Run &run,int Proc_rat ,int& originalcounts  ,  int runID, std::string Eventtype, double cut){
    stringstream ss;
    ss << runID;
    std::string RUN_NUM =  ss.str();
    std::string filepath;//input filepath
    if(Proc_rat == 708)    filepath = "/data/snoplus3/griddata/Prod_7_0_8_9_Preliminary_Scintillator_Gold/ntuples/";
    //if(Proc_rat == 7015)   filepath = "/data/snoplus2/miniPROD_bismsb_newOptics_oldFitter7015/ntuples/";

    else{
        std::cerr << "Error Proc_rat" << std::endl;
        exit(-1);
    }
    std::string input_file  = "ScintFit_2p2"+Eventtype+"Run*"+RUN_NUM+"*.ntuple.root";
    TChain* Tc = new TChain("output"); 

    
    for (const auto &filename : glob( (filepath+input_file).c_str() )) {
        Tc->Add((filename).c_str());
        std::cout<<(filename).c_str()<<std::endl;
    }


    //std::cout<<Tc->GetEntries()<<"   "<<runID<<std::endl;
    
    // check this run number exists or not
    if(Tc->GetEntries() == 0){
        std::cerr << "No such file or Error opening file at run: "<< runID << std::endl;
        delete Tc;

        return;
        //exit(-1);
    }
    
    double mcPosx, mcPosy, mcPosz; int evIndex;
    Tc->SetBranchAddress("mcPosx", &mcPosx);
    Tc->SetBranchAddress("mcPosy", &mcPosy);
    Tc->SetBranchAddress("mcPosz", &mcPosz);
    Tc->SetBranchAddress("evIndex", &evIndex);
    double posx, posy, posz;
    Tc->SetBranchAddress("posx", &posx);
    Tc->SetBranchAddress("posy", &posy);
    Tc->SetBranchAddress("posz", &posz);
    double zoffset = GetZAVoffset(db,run,runID);
    for (int iEntry = 0; iEntry < Tc->GetEntries(); iEntry++){
        Tc->GetEntry(iEntry);
        double mcPoszoffset = mcPosz - zoffset;
        double poszoffset = posz - zoffset;
        double temp_mcR = pow(( mcPosx*mcPosx + mcPosy*mcPosy + mcPoszoffset*mcPoszoffset),0.5 );
        double temp_R = pow(( posx*posx + posy*posy + poszoffset*poszoffset),0.5 );
        if (evIndex <= 0 && temp_mcR < cut){
        //if (evIndex <= 0){
            //mcR.push_back(temp_mcR);
            originalcounts++;
            //std::cout<<originalcounts<<std::endl;
        }
    }
    delete Tc;

    return;
    
}
void GetZFromMySim(int Proc_rat ,int& counts  ,  int runID, std::string Eventtype, double cut){
    stringstream ss;
    ss << runID;
    std::string RUN_NUM =  ss.str();
    std::string filepath;//input filepath
    if(Proc_rat == 708 && Eventtype == "Reactoribd")   filepath = "/data/snoplus3/weiii/antinu/mycuts/Ntuple_sim/Reactoribd/";
    else if(Proc_rat == 708 && Eventtype == "Alphan_Lab_13c")   filepath = "/data/snoplus3/weiii/antinu/mycuts/Ntuple_sim/Alphan_Lab_13c/";
    else if(Proc_rat == 708 && Eventtype == "Geoibd_Th") filepath = "/data/snoplus3/weiii/antinu/mycuts/Ntuple_sim/Geoibd_Th/";
    else if(Proc_rat == 708 && Eventtype == "Geoibd_U") filepath = "/data/snoplus3/weiii/antinu/mycuts/Ntuple_sim/Geoibd_U/";

    else{
        std::cerr << "Error Proc_rat" << std::endl;
        exit(-1);
    }
    std::string input_file  = RUN_NUM+"*.ntuple.root";

    TChain* Tc = new TChain("PromptT"); TChain* DelayTc = new TChain("DelayT"); 
    std::cout<<(filepath+input_file).c_str()<<std::endl;

    for (const auto &filename : glob( (filepath+input_file).c_str() )) {
        Tc->Add(filename.c_str());
        DelayTc->Add(filename.c_str());
    }
    // check this run number exists or not
    if(Tc->GetEntries() == 0 || DelayTc->GetEntries() == 0){
        std::cerr << "No such file or Error opening recon file in run: "<< runID << std::endl;
        //delete Tc; delete DelayTc;
        //exit(-1);
    }
    
    double Prompt_R;double delayedEcorr;//int Prompt_evIndex;
    double dR; double dt;
    Tc->SetBranchAddress("Prompt_R", &Prompt_R);
    //Tc->SetBranchAddress("Prompt_evIndex", &Prompt_evIndex);
    DelayTc->SetBranchAddress("delayedEcorr", &delayedEcorr);
    DelayTc->SetBranchAddress("dR", &dR);
    DelayTc->SetBranchAddress("dt", &dt);
    for (int iEntry = 0; iEntry < Tc->GetEntries(); iEntry++){
        Tc->GetEntry(iEntry);
        double temp_R = Prompt_R; 
        DelayTc->GetEntry(iEntry);
        //if (temp_R < cut && delayedEcorr>1.85 && delayedEcorr<2.4 && dR< 1500 && dt<800000 && dt> 400){
        if (temp_R < cut ){
        //if (temp_R < cut && delayedEcorr>1.85){
            //R.push_back(temp_R);
            counts++;
            //std::cout<<dR<<std::endl;
        }
    }
    delete Tc; delete DelayTc;
    return;
}

std::vector<int>  readrunidfromtable(const std::string& filename){
    std::vector<int> indexes;
    //std::ifstream inputFile("/home/huangp/BiPo/PhysicsRunList/bismsb_bipo214_simlist.txt");
    std::ifstream inputFile(filename);

    if (!inputFile.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return indexes;
    }

    int index;
    while (inputFile >> index) {
        //std::cout << index <<std::endl;
        indexes.push_back(index);
    }

    inputFile.close();
    return indexes;

} 

void effcalculation(int Proc_rat, std::string isotope){
    RAT::DB *db = RAT::DB::Get();
    db->LoadDefaults();
    RAT::DS::Run run;

    std::vector<double> mcR; std::vector<double> R ; std::vector<double> R_4m; std::vector<int> runid;
    int counts = 0; int originalcounts = 0;
    std::string runlist; std::string Eventtype; 
    runlist =  "/home/huangp/AntiNu/preliminary_goldrunlist.txt";
    //runlist =  "/home/huangp/AntiNu/antinu_runlist_UPDATED.txt";
    if(isotope == "Reactoribd"){
        Eventtype = "Reactoribd";
    }
    else if(isotope == "AlphaN"){
        Eventtype = "Alphan_Lab_13c";
    }
    else if(isotope == "Geoibd_U"){
        Eventtype = "Geoibd_U";
    }
    else if(isotope == "Geoibd_Th"){
        Eventtype = "Geoibd_Th";
    }
    else std::cerr << "check isotope argument " << std::endl;
    runid = readrunidfromtable(runlist);
    //std::cout<<"runidsize "<<runid.size()<<std::endl;

    
    for (size_t irun = 0; irun < runid.size(); ++irun){
        GetZFromOriginal(db,run,Proc_rat ,originalcounts  ,  runid[irun], Eventtype, 5700.0);
        //GetZFromMySim(Proc_rat ,counts  ,  runid[irun], Eventtype, 5700.0);
        std::cout<<"Progress: "<<double(irun)/double(runid.size())<<std::endl;
    }
    /*
    std::cout<<"****** Num of Sim in 6m: "<< mcR.size()<< std::endl;
    std::cout<<"****** Num of Recon in 6m: "<< R.size()<< std::endl;
    std::cout<<"****** Eff in 6m: "<<static_cast<float>(R.size())/static_cast<float>(mcR.size())<<std::endl;
    */
    std::cout<<"****** Num of Sim            "<< originalcounts<< std::endl;
    std::cout<<"****** Num of Recon in 5.7m: "<< counts << std::endl;
    std::cout<<"****** Eff in 6m: "<<static_cast<float>(counts)/static_cast<float>(originalcounts)<<std::endl;
    
}


