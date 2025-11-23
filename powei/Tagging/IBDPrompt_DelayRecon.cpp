// compile command:  g++ -g -std=c++1y antinurecon.cpp -o antinurecon.exe `root-config --cflags --libs` -I${RATROOT}/include/libpq -I${RATROOT}/include -I${RATROOT}/include/external -L${RATROOT}/lib -lRATEvent_Linux
// test sim: /home/huangp/AntiNu/antinurecon.exe "/data/snoplus3/griddata/Processing_7_0_8_Preliminary_Scintillator_Gold_300000_308097/ntuples/*300000*.root" 300000  /data/snoplus3/weiii/antinu/mycuts/Ntuple_data/300000.ntuple.root rat-7.0.15 1
// test data: ./antinurecon.exe "/data/snoplus3/griddata/Processing_7_0_8_Preliminary_Scintillator_Gold_300000_308097/ntuples/*306498*.root" 306498 306498.root rat-7.0.8 1

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TRandom.h"
#include "TTree.h"
#include "TChain.h"
#include "TMath.h"
#include <cmath>
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
void printdelayinfo(bool&fitValid, double&posX, double&posY, double&posZ,double delayR, double&posZAfterAVoffset, double&energy,ULong64_t&clock50, ULong64_t&DCapplied,ULong64_t&DCflagged, Int_t&triggerWord, bool is_data);
bool IsDelayEv(bool&fitValid, double&posX, double&posY, double&posZ,double delayR, double&posZAfterAVoffset, double&energy,double&delayedEcorr,ULong64_t&clock50, ULong64_t&DCapplied,ULong64_t&DCflagged, Int_t&triggerWord, bool is_data);
bool IsPromptEv(bool&fitValid, double&posX, double&posY, double&posZ, double promptR, double&posZAfterAVoffset,double&energy, double&promptEcorr, ULong64_t&DCapplied,ULong64_t&DCflagged, Int_t&triggerWord,  bool is_data);
bool Muon_Pileup_Follower_Veto(ULong64_t&DCapplied,ULong64_t&DCflagged,ULong64_t& clock50, int&nhitsCleaned, bool& lone_follower_flag, bool& veto_flag, ULong64_t& lone_start_time, ULong64_t& veto_start_time, ULong64_t& loneFollowerTime, ULong64_t& pileupTime, int& num_vetos);
double EnergyCorrection(const double E, TVector3 pos, const bool is_data, RAT::DU::ReconCalibrator* e_cal);
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
    
    
    if ((DCflagged & 0x80) != 0x80 && nhitsCleaned > 3000){
            //std::cout << "nhitsCleaned: "<< nhitsCleaned <<std::endl;
            //std::cout<< "DCflaggedPass: "<<DCflaggedPass<<std::endl;
            
        
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



double EnergyCorrection(const double E, TVector3 pos, const bool is_data, RAT::DU::ReconCalibrator* e_cal) {
    // Data vs MC energy correction (Tony's)
    double Ecorr = e_cal->CalibrateEnergyRTF(is_data, E, std::sqrt(pos.X()*pos.X() + pos.Y()*pos.Y()), pos.Z()); // gives the new E

   

    return Ecorr;
}





bool IsDelayEv(bool&fitValid, double&posX, double&posY, double&posZ,double delayR, double&posZAfterAVoffset, double&energy,double&delayedEcorr,ULong64_t&clock50, ULong64_t&DCapplied,ULong64_t&DCflagged, Int_t&triggerWord, bool is_data){
    // 0 : not pass all the cuts; 1: is a delay candidates

    // orphan check
    if (is_data == 1){
        if (triggerWord == 0) return false;
    }
    if(!fitValid ) return false;
    
    // FV check    
    if( delayR > 5700.0 ) return false;
    //std::cout<< "R: "<<R <<std::endl;
    // apply DC mask
    if (is_data == 1){
        if (((DCapplied & 0x2100000042C2) & DCflagged ) != (DCapplied & 0x2100000042C2)) return false;
    }
    //std::cout<< "pass the DC"<<std::endl;
    // apply delayed energy cuts
    if(delayedEcorr < 1.7 || delayedEcorr > 2.5)  return false;
    //if(energy < 1.85 || energy > 2.5)  return false;
    //std::cout<< "pass the energy"<<std::endl;
    
    return true;
}
void printdelayinfo(bool&fitValid, double&posX, double&posY, double&posZ,double delayR, double&posZAfterAVoffset, double&energy,ULong64_t&clock50, ULong64_t&DCapplied,ULong64_t&DCflagged, Int_t&triggerWord, bool is_data){
    std::cout <<"fitValid "<< fitValid<<std::endl;
    std::cout <<"delayR "<< delayR<<std::endl;
    std::cout <<" posZAfterAVoffset "<< posZAfterAVoffset<<std::endl;
    std::cout <<" energy "<< energy<<std::endl;
    std::cout <<" triggerWord "<<triggerWord <<std::endl;
    
}
bool IsPromptEv(bool&fitValid, double&posX, double&posY, double&posZ, double promptR, double&posZAfterAVoffset,double&energy, double&promptEcorr, ULong64_t&DCapplied,ULong64_t&DCflagged, Int_t&triggerWord,  bool is_data){

    // orphan check
    if (is_data == 1){
        if (triggerWord == 0) return false;
    }
    // fit valid check
    //std::cout<<fitValid<<std::endl;

    if(!fitValid ) return false; 
    
    // dR check
    //std::cout<<dR<<endl;

    //std::cout<<"pass DR cuts"<<std::endl;
    if (is_data == 1){
        if (((DCapplied & 0x2100000042C2) & DCflagged ) != (DCapplied & 0x2100000042C2)) return false;
    }
    if( promptR> 5700.0 ) return false;
    //std::cout<< "pass the prompt radial cuts"<<std::endl;
    //std::cout<< "energy"<< energy<<std::endl;
    if(promptEcorr < 0.9 || promptEcorr > 8.0) return false;
    //if(energy < 0.9 || energy > 8.0) return false;
    //std::cout<< "pass the energy cuts"<<std::endl;
    //std::cout<< "pass promptcuts"<<std::endl;
    return true;
    
}

    


void IBDACCidentalRecon(std::string inputpath, int RUN_NUM, std::string output_root_address, std::string Proc_rat, bool is_data){
    TChain *Tc = new TChain("output");
    TChain *PoTc = new TChain("PoT"); TChain *BiTc = new TChain("BiT");
    
    std::string filepath;//input filepath
    if(is_data == 0 && Proc_rat == "rat-7.0.8"){
        //filepath = "/data/snoplus3/griddata/Prod_7_0_8_9_Preliminary_Scintillator_Gold/ntuples";
        //std::string input_file  = "*"+datatype+"*.ntuple.root";
        
        Tc->Add(inputpath.c_str());
    }
    else if(is_data == 1 && Proc_rat == "rat-7.0.8"){
        BiTc->Add(("/data/snoplus3/weiii/BiPo214/fullFill/rat-7.0.8/Ntuple/*"+std::to_string(RUN_NUM)+"*").c_str());
        PoTc->Add(("/data/snoplus3/weiii/BiPo214/fullFill/rat-7.0.8/Ntuple*"+std::to_string(RUN_NUM)+"*").c_str());
        Tc->Add(inputpath.c_str());
        //filepath = inputpath;
        //std::string input_file  = "Analysis*_r0000"+std::to_string(RUN_NUM)+"*.ntuple.root";

        //std::string filename = filepath + "/"+input_file;

        //Tc->Add(filename.c_str());
        
        // check this run number exists or not
        if(Tc->GetEntries() == 0){
            std::cerr << "Error opening files or files corruption"  << std::endl;
            exit(-1);
        }

        if(BiTc->GetEntries() == 0 || PoTc->GetEntries() == 0){
            std::cerr << "Error opening BiPo tagging files or corruption" << std::endl;
            //exit(-1);
        }
    }
    else{
        std::cerr << "Error Processing RAT Argument" << std::endl;
        exit(-1);
    }
    std::cout<< "Ready to process files in "<< filepath <<std::endl; 
    
     // define BiPo variables to load in
    int Po_eventid; int Bi_eventid;
    BiTc->SetBranchAddress("Bi_eventid", &Bi_eventid);
    PoTc->SetBranchAddress("Po_eventid", &Po_eventid);
    std::vector<int> eventid;
    for(int lEntry = 0; lEntry< BiTc->GetEntries(); lEntry++){
        BiTc->GetEntry(lEntry); PoTc->GetEntry(lEntry);
        eventid.push_back(Bi_eventid);
        eventid.push_back(Po_eventid);
    }
    std::sort(eventid.begin(), eventid.end());

// define variables to save the event information
    double posX, posY, posZ, energy, delayedEcorr, promptEcorr,delayposZAfterAVoffset, promptposZAfterAVoffset;
    int nhits; int nhitsCleaned; Double_t correctedNhits;
    bool scintFit, fitValid, partialFit;
    ULong64_t clock50;
    ULong64_t DCapplied;
    ULong64_t DCflagged;
    Int_t triggerWord;
    Int_t owlnhits;
    double promptR; double delayR;
    int ncandidates; // num of Bi per Po
    int gtid; int runID;

    //initialize all parameters to trck livetime
    bool lone_follower_flag = false; bool veto_flag = false; 
    ULong64_t lone_start_time=0.; ULong64_t veto_start_time=0.; 
    ULong64_t loneFollowerTime=0.; ULong64_t pileupTime=0.;
    int num_vetos=0;int num_spallationvetos=0;
    int64_t delayedTime, promptTime;
    double highNhitDelay, owlNhitDelay;
    int64_t highNhitTime = -99999999999;
    int64_t owlNhitTime = -99999999999;
    // Create a output file, TTree and branches
    //std::string output_root_address; //output filepath
    
    
   


    TFile *f_out = new TFile( output_root_address.c_str(),"RECREATE");
    TTree *DelayT = new TTree("DelayT","Po After Cuts TTree");
    TTree *PromptT = new TTree("PromptT","Bi After Cuts TTree");
    

    DelayT->Branch("Delay_fitValid", &fitValid);
    DelayT->Branch("Delay_scintFit", &scintFit);
    DelayT->Branch("Delay_partialFit", &partialFit);
    DelayT->Branch("Delay_posX", &posX);
    DelayT->Branch("Delay_posY", &posY);
    DelayT->Branch("Delay_posZ", &posZ);
    DelayT->Branch("Delay_posZAfterAVoffset", &delayposZAfterAVoffset);
    DelayT->Branch("Delay_clockCount50", &clock50);
    DelayT->Branch("Delay_energy", &energy);
    DelayT->Branch("Delay_R", &delayR);
    DelayT->Branch("Delay_nhits", &nhits);
    DelayT->Branch("Delay_nhitsCleaned", &nhitsCleaned);
    DelayT->Branch("Delay_correctedNhits", &correctedNhits);
    DelayT->Branch("ncandidates", &ncandidates);
    DelayT->Branch("Delay_eventid", &gtid);
    DelayT->Branch("Delay_triggerWord", &triggerWord);
    DelayT->Branch("Delay_owlnhits", &owlnhits);
    DelayT->Branch("Delay_runID", &runID);
    DelayT->Branch("delayedEcorr",&delayedEcorr);

    

    PromptT->Branch("Prompt_fitValid", &fitValid);
    PromptT->Branch("Prompt_scintFit", &scintFit);
    PromptT->Branch("Prompt_partialFit", &partialFit);
    PromptT->Branch("Prompt_posX", &posX);
    PromptT->Branch("Prompt_posY", &posY);
    PromptT->Branch("Prompt_posZ", &posZ);
    PromptT->Branch("Prompt_posZAfterAVoffset", &promptposZAfterAVoffset);
    PromptT->Branch("Prompt_clockCount50", &clock50);
    PromptT->Branch("Prompt_energy", &energy);
    PromptT->Branch("Prompt_R", &promptR);
    PromptT->Branch("Prompt_nhits", &nhits);
    PromptT->Branch("Prompt_nhitsCleaned", &nhitsCleaned);
    PromptT->Branch("Prompt_correctedNhits", &correctedNhits);
    PromptT->Branch("Prompt_eventid", &gtid);
    PromptT->Branch("Prompt_triggerWord", &triggerWord);
    PromptT->Branch("Prompt_owlnhits", &owlnhits);
    PromptT->Branch("Prompt_runID", &runID);
    PromptT->Branch("promptEcorr",&promptEcorr);
    




    // to check the reconstruction was a) run at all and b) converged to a good solution
    Tc->SetBranchAddress("scintFit", &scintFit);
    Tc->SetBranchAddress("fitValid", &fitValid);
    Tc->SetBranchAddress("partialFit", &partialFit);
    // get the reconstructed position to check event falls within FV and work out the dR between Bi and Po candidate
    Tc->SetBranchAddress("posx", &posX);
    Tc->SetBranchAddress("posy", &posY);
    Tc->SetBranchAddress("posz", &posZ);

    // We use the clockCount50 (50 MHz) clock to find the inter-event time
    Tc->SetBranchAddress("clockCount50", &clock50);

    // Use the reconstructed energy to tag the Po and then the Bi
    Tc->SetBranchAddress("energy", &energy);
    Tc->SetBranchAddress("nhits", &nhits);
    Tc->SetBranchAddress("owlnhits", &owlnhits);
    Tc->SetBranchAddress("nhitsCleaned", &nhitsCleaned);
    Tc->SetBranchAddress("correctedNhits", &correctedNhits);

    // DC mask variables
    Tc->SetBranchAddress("dcApplied", &DCapplied);
    Tc->SetBranchAddress("dcFlagged", &DCflagged);
    Tc->SetBranchAddress("eventID", &gtid);
    Tc->SetBranchAddress("runID", &runID);

    // trigger word
    Tc->SetBranchAddress("triggerWord", &triggerWord);

    int Pocounts = 0;
    double ZAVoffset;
    // get AV offset in zaxis to do correct R=x^2+y^2+z^2 calculation
    ZAVoffset = GetZAVoffset(RUN_NUM);


    // induce energy correction
    RAT::DU::ReconCalibrator* e_cal = RAT::DU::ReconCalibrator::Get();
    
    // now iterate through every event in the ntuple to look for the DELAYED (Po) events
    for (int iEntry = 0; iEntry < Tc->GetEntries(); iEntry++){
        
        
        ncandidates = 0; // refresh the Bi count
        Tc->GetEntry(iEntry);
        // apply muon Veto
        if (is_data == 1){
            if (nhits > 3000) highNhitTime = delayedTime;
            highNhitDelay = ((delayedTime - highNhitTime) & 0x7FFFFFFFFFF) / 50E6;  // [s] dealing with clock rollover

            if( !Muon_Pileup_Follower_Veto(DCapplied,DCflagged,clock50, nhitsCleaned, lone_follower_flag, veto_flag, lone_start_time, veto_start_time, loneFollowerTime, pileupTime, num_vetos)){
                continue;
            }
        }
        delayedTime = int64_t(clock50);
        if (owlnhits > 3) owlNhitTime = delayedTime;
        owlNhitDelay = ((delayedTime - owlNhitTime) & 0x7FFFFFFFFFF) / 50.0;  // [us] dealing with clock rollover
        if (owlNhitDelay < 10){
            continue;
        }

        //BiPo veto
        if (is_data == 1){
            if ( std::find(eventid.begin(), eventid.end(), gtid) != eventid.end() ) continue;

        }
        

        delayposZAfterAVoffset = posZ - ZAVoffset;
        delayR = pow(( posX*posX + posY*posY + delayposZAfterAVoffset*delayposZAfterAVoffset),0.5 );
        
        TVector3 delayedPos = TVector3(posX, posY, posZ);
        delayedEcorr = EnergyCorrection(energy, delayedPos, is_data, e_cal);
        if (!IsDelayEv(fitValid, posX, posY, posZ, delayR,delayposZAfterAVoffset, energy,delayedEcorr, clock50,DCapplied,DCflagged,triggerWord, is_data)){
            continue;
        }
        
        ULong64_t Delay_clock50 = clock50;
        
        
        DelayT->Fill();
    }

    for (int jEntry = 0; jEntry < Tc->GetEntries(); jEntry++){
        Tc->GetEntry(jEntry);
        promptTime = int64_t(clock50);
        if (is_data == 1){
            highNhitDelay = ((promptTime - highNhitTime) & 0x7FFFFFFFFFF) / 50E6;  // [s] dealing with clock rollover
            if (highNhitDelay < 20) break;
        
        }

        
        if (owlnhits > 3) owlNhitTime = promptTime;
        owlNhitDelay = ((promptTime - owlNhitTime) & 0x7FFFFFFFFFF) / 50.0;  // [us] dealing with clock rollover
        if (owlNhitDelay < 10){
            continue;
        }
        
        //BiPo Veto
        if (is_data == 1){
            if ( std::find(eventid.begin(), eventid.end(), gtid) != eventid.end() ) continue;
        }

        promptposZAfterAVoffset = posZ - ZAVoffset;
        promptR = pow(( posX*posX + posY*posY + promptposZAfterAVoffset*promptposZAfterAVoffset),0.5 );

        //std::cout<<jEntry<<" IsPromptEv: "<<IsPromptEv(fitValid, posX, posY, posZ, posZAfterAVoffset,energy,dR, DCapplied,DCflagged, triggerWord)<<std::endl;
        TVector3 promptPos = TVector3(posX, posY, posZ);
            promptEcorr = EnergyCorrection(energy, promptPos, is_data, e_cal);
        if (!IsPromptEv(fitValid, posX, posY, posZ, promptR,promptposZAfterAVoffset,energy, promptEcorr,DCapplied,DCflagged, triggerWord, is_data)){
            continue;
        }
        
        PromptT->Fill();
        
        
        

        
    } 
        
        
        

    
    f_out->Write();
    f_out->Close();
    //std::cout<<"Entries:  "<<Tc->GetEntries()<<std::endl;

    std::cout<<"End of Coincidence tagging"<<std::endl;


    
    
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
    
    IBDACCidentalRecon(inputfilepath, RUN_NUM,output_root_address,ratversion,is_data);


    return 0;
}
