
#include<stdlib.h>
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
std::vector<std::string> glob(const char *);
double GetZAVoffset(int);
bool Muon_Pileup_Follower_Veto(ULong64_t&DCapplied,ULong64_t&DCflagged,ULong64_t& clock50, int&nhitsCleaned, bool& lone_follower_flag, bool& veto_flag, ULong64_t& lone_start_time, ULong64_t& veto_start_time, ULong64_t& loneFollowerTime, ULong64_t& pileupTime, int& num_vetos);
bool IsDelayEv(int run, bool&fitValid, double&posX, double&posY, double&posZ,double delayR, double&posZAfterAVoffset, Double_t& correctedNhits,double&energy,ULong64_t&clock50, ULong64_t&DCapplied,ULong64_t&DCflagged, Int_t&triggerWord, bool is_data);
bool IsPromptEv(int run,bool&fitValid, double&posX, double&posY, double&posZ, double promptR, double&posZAfterAVoffset,Double_t& correctedNhits,double&energy,double&dR, ULong64_t&DCapplied,ULong64_t&DCflagged, Int_t&triggerWord,  bool is_data);
// use glob to match those files passing the wildcard
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

double GetZAVoffset(int runID){
    RAT::DB *db = RAT::DB::Get();
    db->LoadDefaults();
    RAT::DS::Run run;
    run.SetRunID(runID);
    db->BeginOfRun(run);
    std::vector<double> AVPos = RAT::GeoUtil::UpdateAVOffsetVectorFromDB();
    double zOff = AVPos[2];
    std::cout << "AV Shift is: " << zOff << std::endl;
    return zOff;
}




bool Muon_Pileup_Follower_Veto(ULong64_t&DCapplied,ULong64_t&DCflagged,ULong64_t& clock50, int&nhitsCleaned, bool& lone_follower_flag, bool& veto_flag, ULong64_t& lone_start_time, ULong64_t& veto_start_time, ULong64_t& loneFollowerTime, ULong64_t& pileupTime, int& num_vetos){
    if ((DCflagged & 0x80) != 0x80 && nhitsCleaned > 5000){
            //std::cout << "nhitsCleaned: "<< nhitsCleaned <<std::endl;
            
        
          // check if we are inside a lone muon follower veto window
          if (lone_follower_flag == true){
              // add the dT and switch off the muon follower veto
              ULong64_t deltaT = ((clock50-lone_start_time) & 0x7FFFFFFFFFF)*20.0;
              loneFollowerTime += deltaT;
              lone_follower_flag = false;
          }

          if (veto_flag == false){
              num_vetos++;
              veto_flag = true;
              veto_start_time = clock50;
              return false;
          }
          else{
              // we have pileup! Need to calculate the additional pileup time.
              
              ULong64_t deltaT = ((clock50-veto_start_time) & 0x7FFFFFFFFFF)*20.0;
              std::cout<<"we have pielup!!! "<<"increment pileuptime: "<<deltaT<<std::endl;
              pileupTime += deltaT;

              // reset the veto window
              veto_start_time = clock50;
              return false;
          }
      }

      // now handle the veto window for follower events from highE / muon veto
      if (veto_flag == true){
          ULong64_t deltaT = ((clock50-veto_start_time) & 0x7FFFFFFFFFF)*20.0;

          // check if event falls inside veto window
          if (deltaT < 2.0e10){
            return false;
          }
          else{
            // we are no longer inside the veto time window --> switch it off!
            veto_flag = false;
          }
      }

      // check if we have a lone muon follower (i.e. muon at end of previous run)
      // this can only trigger if there was a muon at the end of the previous run
      if ((DCflagged&0x4000)!=0x4000){
          if (lone_follower_flag == false){
            lone_follower_flag = true; 
            lone_start_time = clock50;
            return false;
          } else{
            // we are within a lone follower veto period!
            ULong64_t deltaT = ((clock50 - lone_start_time) & 0x7FFFFFFFFFF)*20.0;
            loneFollowerTime += deltaT;
            lone_start_time = clock50;
            return false;
          }
      }
      else{
        if(lone_follower_flag == true){
            lone_follower_flag = false;
        }
      }
    return true;
}
bool IsDelayEv(int run, bool&fitValid, double&posX, double&posY, double&posZ,double delayR, double&posZAfterAVoffset, Double_t& correctedNhits,double&energy,ULong64_t&clock50, ULong64_t&DCapplied,ULong64_t&DCflagged, Int_t&triggerWord, bool is_data){
    // 0 : not pass all the cuts; 1: is a delay candidates

    // orphan check
    if (is_data == 1){
        if (triggerWord == 0) return false;
    }
    if(!fitValid ) return false;
    
    // FV check    
    if( delayR > 6000.0 ) return false;
    //std::cout<< "R: "<<R <<std::endl;
    // apply DC mask
    if (is_data == 1){
        if (((DCapplied & 0x2100000042C2) & DCflagged ) != (DCapplied & 0x2100000042C2)) return false;
    }
    //std::cout<< "pass the DC"<<std::endl;
    // apply delayed energy cuts
    if(run>=350000){
            //if(energy < 1.36 || energy > 2.1)  return false;
            
            //if(correctedNhits < 425 || correctedNhits > 530)   return false;
        }
    else{
            if(energy < 0.6 || energy > 2.5)  return false;
    }
    //std::cout<< "Po pass the Energy"<<std::endl;
    return true;
}
bool IsPromptEv(int run,bool&fitValid, double&posX, double&posY, double&posZ, double promptR, double&posZAfterAVoffset,Double_t& correctedNhits,double&energy,double&dR, ULong64_t&DCapplied,ULong64_t&DCflagged, Int_t&triggerWord,  bool is_data){

    // orphan check
    if (is_data == 1){
        if (triggerWord == 0) return false;
    }

    if(!fitValid ) return false; 
    // dR check
    if(dR > 800) return false; 
   // std::cout<< "pass the dR cuts"<<std::endl;
    // pass DC mask
    //std::cout<<DCapplied<<endl;
    if (is_data == 1){
        if (((DCapplied & 0x2100000042C2) & DCflagged ) != (DCapplied & 0x2100000042C2)) return false;
    }
    if( promptR > 6000.0 ) return false;
    //std::cout<< "pass the prompt radial cuts"<<std::endl;
    if(run>=350000){
                //std::cout<< "Energy: "<<energy<<std::endl;
        //if(energy < 1.4 || energy > 4.0)  return false;
        //if(correctedNhits < 350 || correctedNhits > 850)  return false;
    }
    else{
        if(energy < 0.58 || energy > 1.08) return false;
    }
    //std::cout<< "Bi pass the Energy"<<std::endl;
    return true;
    
}
void RnPo215Recon(std::string inputpath, int RUN_NUM, std::string output_root_address, std::string Proc_rat, bool is_data){// detector_status: 0 is water phase, 1 is partial fill, 2 is  fulfill
    TChain *Tc = new TChain("output");
    if(is_data == 0 && Proc_rat == "rat-7.0.15") Tc->Add(inputpath.c_str());
    else if(is_data == 1 && Proc_rat == "rat-7.0.8")    Tc->Add(inputpath.c_str());//filepath = "/data/snoplus3/SNOplusData/processing/fullFill/rat-7.0.8/ntuples/";
    else if(is_data == 1 && Proc_rat == "rat-7.0.15")   Tc->Add(inputpath.c_str());//filepath = "/data/snoplus3/SNOplusData/processing/fullFill/rat-7.0.15/ntuples/";

    else{
        std::cerr << "Error InputFile" << std::endl;
        exit(-1);
    }
    //std::string input_file  = "Analysis*_r0000"+RUN_NUM+"*.ntuple.root";

    
    // check this run number exists or not
    if(Tc->GetEntries() == 0){
        std::cerr << "Error opening file" << std::endl;
        exit(-1);
    }
    std::cout<< "Ready to process files in "<< inputpath <<std::endl; 


    // define variables to save the event information
    double posX, posY, posZ, energy, delayposZAfterAVoffset, promptposZAfterAVoffset;
    int nhits; int nhitsCleaned; Double_t correctedNhits;
    bool scintFit, fitValid, partialFit;
    ULong64_t clock50;
    ULong64_t DCapplied;
    ULong64_t DCflagged;
    Int_t triggerWord;
    double promptR; double delayR;
    double dR; // record distance difference between Po/Bi
    double dt;// time difference between Po/Bi
    int ncandidates; // num of Bi per Po
    int gtid;
    int evIndex;
    int pdg1; int pdg2;
    //initialize all parameters to trck livetime
    bool lone_follower_flag = false; bool veto_flag = false; 
    ULong64_t lone_start_time=0.; ULong64_t veto_start_time=0.; 
    ULong64_t loneFollowerTime=0.; ULong64_t pileupTime=0.;
    int num_vetos=0; 
    // Create a output file, TTree and branches
    std::string outputfilepath; //output filepath
    TFile *f_out = new TFile( output_root_address.c_str(),"RECREATE");
    
    TTree *PoT = new TTree("PoT","Po After Cuts TTree");
    TTree *BiT = new TTree("BiT","Bi After Cuts TTree");

    
    

    PoT->Branch("Po_fitValid", &fitValid);
    PoT->Branch("Po_scintFit", &scintFit);
    PoT->Branch("Po_partialFit", &partialFit);
    PoT->Branch("Po_posX", &posX);
    PoT->Branch("Po_posY", &posY);
    PoT->Branch("Po_posZ", &posZ);
    PoT->Branch("Po_posZAfterAVoffset", &delayposZAfterAVoffset);
    PoT->Branch("Po_clockCount50", &clock50);
    PoT->Branch("Po_energy", &energy);
    PoT->Branch("Po_R", &delayR);
    PoT->Branch("Po_nhits", &nhits);
    PoT->Branch("Po_correctedNhits", &correctedNhits);
    PoT->Branch("Po_nhitsCleaned", &nhitsCleaned);
    PoT->Branch("ncandidates", &ncandidates);
    PoT->Branch("Po_eventid", &gtid);
    PoT->Branch("Po_triggerWord", &triggerWord);
    PoT->Branch("Po_evIndex", &evIndex);
    PoT->Branch("Po_pdg1", &pdg1);
    PoT->Branch("Po_pdg2", &pdg2);
    
    

    BiT->Branch("Bi_fitValid", &fitValid);
    BiT->Branch("Bi_scintFit", &scintFit);
    BiT->Branch("Bi_partialFit", &partialFit);
    BiT->Branch("Bi_posX", &posX);
    BiT->Branch("Bi_posY", &posY);
    BiT->Branch("Bi_posZ", &posZ);
    BiT->Branch("Bi_posZAfterAVoffset", &promptposZAfterAVoffset);
    BiT->Branch("Bi_clockCount50", &clock50);
    BiT->Branch("Bi_energy", &energy);
    BiT->Branch("dR", &dR);
    BiT->Branch("dt", &dt);
    BiT->Branch("Bi_R", &promptR);
    BiT->Branch("Bi_nhits", &nhits);
    BiT->Branch("Bi_nhitsCleaned", &nhitsCleaned);
    BiT->Branch("Bi_correctedNhits", &correctedNhits);
    BiT->Branch("Bi_eventid", &gtid);
    BiT->Branch("Bi_triggerWord", &triggerWord);
    BiT->Branch("Bi_evIndex", &evIndex);
    




    // to check the reconstruction was a) run at all and b) converged to a good solution
    Tc->SetBranchAddress("scintFit", &scintFit);
    Tc->SetBranchAddress("fitValid", &fitValid);
    Tc->SetBranchAddress("partialFit", &partialFit);
    // get the reconstructed position to check event falls within FV and work out the dR between Bi and Po candidate
    Tc->SetBranchAddress("posx", &posX);
    Tc->SetBranchAddress("posy", &posY);
    Tc->SetBranchAddress("posz", &posZ);

    // We use the clockCount50 (50 MHz) clock to find the inter-event tim...posZAfterAVoffset￼✕Bi_clockCount50￼✕

    Tc->SetBranchAddress("clockCount50", &clock50);

    // Use the reconstructed energy to tag the Po and then the Bi
    Tc->SetBranchAddress("energy", &energy);
    Tc->SetBranchAddress("nhits", &nhits);
    Tc->SetBranchAddress("nhitsCleaned", &nhitsCleaned);
    Tc->SetBranchAddress("correctedNhits", &correctedNhits);


    // DC mask variables
    Tc->SetBranchAddress("dcApplied", &DCapplied);
    Tc->SetBranchAddress("dcFlagged", &DCflagged);
    Tc->SetBranchAddress("eventID", &gtid);
    Tc->SetBranchAddress("evIndex", &evIndex);
    Tc->SetBranchAddress("pdg1", &pdg1);
    Tc->SetBranchAddress("pdg2", &pdg2);

    // trigger word
    Tc->SetBranchAddress("triggerWord", &triggerWord);


    int Pocounts = 0;
    int orphancounts = 0;
    double ZAVoffset;
    // get AV offset in zaxis to do correct R=x^2+y^2+z^2 calculation
    ZAVoffset = GetZAVoffset(RUN_NUM);

    // now iterate through every event in the ntuple to look for the DELAYED (Po) events
    for (int iEntry = 0; iEntry < Tc->GetEntries(); iEntry++){
        ncandidates = 0; // refresh the Bi count
        Tc->GetEntry(iEntry);
        
        // apply muon Veto
        if (is_data == 1){
            if( !Muon_Pileup_Follower_Veto(DCapplied,DCflagged,clock50, nhitsCleaned, lone_follower_flag, veto_flag, lone_start_time, veto_start_time, loneFollowerTime, pileupTime, num_vetos)) continue;
        }
        delayposZAfterAVoffset = posZ - ZAVoffset;
        delayR = pow(( posX*posX + posY*posY + delayposZAfterAVoffset*delayposZAfterAVoffset),0.5 );

        //std::cout<<"IsDelayEv "<<IsDelayEv(fitValid, posX, posY, posZ, posZAfterAVoffset, energy, clock50,DCapplied,DCflagged,triggerWord)<<std::endl;
        if (!IsDelayEv(RUN_NUM,fitValid, posX, posY, posZ, delayR,delayposZAfterAVoffset, correctedNhits,energy, clock50,DCapplied,DCflagged,triggerWord, is_data)){
            continue;
        }
        
        
        
        
        ULong64_t Delay_clock50 = clock50;

        
        
        // record x y z of Po
        double delayX = posX; double delayY =  posY ; double delayZ = delayposZAfterAVoffset;
        
        // if successful, loop backwards in time to find the Bi candidate
        for (int jEntry = iEntry - 1; jEntry >= 0; jEntry --){
            

            Tc->GetEntry(jEntry);
            ULong64_t Prompt_clock50 = clock50;
            dt = ((Delay_clock50 - Prompt_clock50) & 0x7FFFFFFFFFF)* 20.0/1000.0;//us
            //std::cout<<"dt: "<<dt<<std::endl;
            if( dt < 1100) continue;
            if( dt > 11000 ) break;
            
    
            // apply dR in mm
            promptposZAfterAVoffset = posZ - ZAVoffset;

            dR = pow( (posX- delayX)*(posX- delayX)+ (posY- delayY)*(posY- delayY)+ (promptposZAfterAVoffset- delayZ)*(promptposZAfterAVoffset- delayZ) ,0.5);
            //std::cout<< "dR: "<<dR<<std::endl;

            promptR = pow(( posX*posX + posY*posY + promptposZAfterAVoffset*promptposZAfterAVoffset),0.5 );


            //std::cout<<jEntry<<" IsPromptEv: "<<IsPromptEv(fitValid, posX, posY, posZ, posZAfterAVoffset,energy,dR, DCapplied,DCflagged, triggerWord)<<std::endl;
            if (!IsPromptEv(RUN_NUM,fitValid, posX, posY, posZ, promptR,promptposZAfterAVoffset,correctedNhits,energy,dR, DCapplied,DCflagged, triggerWord, is_data)){
                continue;
            }

            ncandidates++;
            
            
            BiT->Fill();
            Tc->GetEntry(iEntry);
            PoT->Fill();
            
            
        } 
        
        

    }
    f_out->Write();  
    
    // calculate lifetime and save results into file
    

    
}

int main(int argc, char** argv) {
    std::string inputfilepath = argv[1];
    int RUN_NUM = std::stoi(argv[2]);
    std::string output_root_address = argv[3];
    std::string ratversion = argv[4];  // Old legacy entry, left for backwards compatibility with wrapper code (unused here)
    //int runid = std::stoi(argv[3]);
    bool is_data = std::stoi(argv[5]); 
    
    
    // Addresses of simulation output files to be analysed
    std::vector<std::string> input_files;
    RnPo215Recon(inputfilepath, RUN_NUM,output_root_address,ratversion,is_data);


    return 0;
}
