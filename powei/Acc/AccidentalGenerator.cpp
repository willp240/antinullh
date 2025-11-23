// compile command:  g++ -g -std=c++1y AccidentalGenerator.cpp -o AccidentalGenerator.exe `root-config --cflags --libs` -I${RATROOT}/include/libpq -I${RATROOT}/include     -I${RATROOT}/include/external -L${RATROOT}/lib -lRATEvent_Linux             
// test data: ./AccidentalGenerator.exe 2p2goldlist.txt test.root rat-7.0.8 1000
#include "TFile.h"
#include "TNtuple.h"
#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TRandom3.h"
#include "TTree.h"
#include "TChain.h"
#include "TMath.h"
#include <time.h>

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <glob.h>
#include <RAT/DU/DSReader.hh>                                               
#include <RAT/DU/Utility.hh>
#include <RAT/DS/Entry.hh>
#include <RAT/GeoUtils.hh>

// generate a specifified number of random variables
std::vector<double> randomnumber_generator(int N){
    
    std::vector<double> randomarr;
    TRandom *r3 = new TRandom3(time(NULL));
    for (int i=0;i<N;i++) {
    //std::cout<<r3->Rndm(i)<<std::endl;
    randomarr.push_back(r3->Rndm(i));
    }
    return randomarr;
}
std::vector<std::string> loadrunidfromlist(std::string list_name){
    std::vector<std::string> runidlist;
    std::ifstream inputFile(list_name); 
    if (!inputFile) {
        std::cerr << "Unable to open file" << std::endl;
        exit(-1);
    }

    std::string line;
    while (std::getline(inputFile, line)) {
        runidlist.push_back(line);
    }

    inputFile.close();
    return runidlist;
}

std::vector<std::string> loadsubrunidfromprocesslist(std::string list_name){
    std::vector<std::string> runidlist;
    std::ifstream inputFile(list_name); 
    if (!inputFile) {
        std::cerr << "Unable to open file" << std::endl;
        exit(-1);
    }

    std::string line;
    while (std::getline(inputFile, line)) {
        //std::cout<<line.substr(0,45)<<std::endl;

        runidlist.push_back(line.substr(0,45));
    }

    inputFile.close();
    return runidlist;
}


void generate_accidental( std::string runlist,std::string output_root_address, std::string rat_ver, int N){
    // generate random numbers
    std::vector<double> randomarr = randomnumber_generator(4*N);
    // set the output file, tree
    TFile *f_out = new TFile( output_root_address.c_str(),"RECREATE");
    TTree *AccT = new TTree("AccT","Tree of Accidental");
    

    
    // load prompt and delay root files
    std::string PATH = "/data/snoplus3/weiii/antinu/mycuts/Ntuple_data/accidental/";
    //std::vector<std::string> runidlist = loadrunidfromlist(runlist);
    std::vector<std::string> subrunidlist = loadsubrunidfromprocesslist(runlist);
    
    int subfilesize = subrunidlist.size();

    // set the variables we want to save (SetBranchAddress, Branch)
    int prompt_gtid; double prompt_energy; double promptEcorr;ULong64_t prompt_clock50; double prompt_posx; double prompt_posy; double prompt_posz ; double promptR; 
    int delay_gtid ;double delay_energy; double delayedEcorr;ULong64_t delay_clock50; double delay_posx; double delay_posy; double delay_posz ;double delayR;
    double dR; double dt; double dt_model; int Delay_runID;
    

    AccT->Branch("dR",&dR);
    AccT->Branch("dt",&dt);
    AccT->Branch("dt_model",&dt_model);
    AccT->Branch("promptR",&promptR);
    AccT->Branch("delayR",&delayR);
    AccT->Branch("promptE",&prompt_energy);
    AccT->Branch("promptEcorr",&promptEcorr);
    AccT->Branch("delayE",&delay_energy);
    AccT->Branch("delayedEcorr",&delayedEcorr);
    AccT->Branch("prompt_gtid",&prompt_gtid);
    AccT->Branch("delay_gtid",&delay_gtid);
    AccT->Branch("RunID",&Delay_runID);

    for( int i = 0; i < randomarr.size(); i+=4){

        // load prompt and delay root files
        std::string input_file = subrunidlist[(int)round(randomarr[i]*subfilesize)];
        TFile* inputFile = TFile::Open((PATH+"AfterCuts_"+input_file).c_str(), "READ");
        if (!inputFile || inputFile->IsZombie()) {
            std::cerr << "Error: could not open input file!" << std::endl;
            continue;
        }

        TTree* DelayTc = (TTree*)inputFile->Get("DelayT");
        TTree* PromptTc = (TTree*)inputFile->Get("PromptT");
        if (!DelayTc || !PromptTc) {
            std::cerr << "Error: TTree DelayT or PromptT not found!" << std::endl;
            inputFile->Close();
            continue;
        }
        /*
        TChain *PromptTc;TChain *DelayTc;
        //PromptTc->Add((PATH+"*"+input_file).c_str()); DelayTc->Add((PATH+"*"+input_file).c_str());
        
        //if(PromptTc->GetEntries() == 0 || DelayTc->GetEntries() == 0) continue;
        try {
            std::string input_file = subrunidlist[(int)round(randomarr[i]*subfilesize)];
            PromptTc = new TChain("PromptT"); DelayTc = new TChain("DelayT");

            //std::cout<<PATH+"*"+input_file<<std::endl;
            PromptTc->Add((PATH+"*"+input_file).c_str()); DelayTc->Add((PATH+"*"+input_file).c_str());
            if(PromptTc->GetEntries() == 0 || DelayTc->GetEntries() == 0) throw(input_file);
        }
        catch(std::string input_file){
            std::cout<<PATH+"*"+input_file+ "cannot be found or no entry in it"<<std::endl;
            std::cout<<PromptTc->GetEntries() <<" "<< DelayTc->GetEntries()<<std::endl;
            continue;
        }
        
        */
        //std::cout<<"found a file"<<std::endl;
        
        PromptTc->SetBranchAddress("promptEcorr",&promptEcorr); 
        PromptTc->SetBranchAddress("Prompt_energy",&prompt_energy); 
        PromptTc->SetBranchAddress("Prompt_clockCount50",&prompt_clock50); ; 
        PromptTc->SetBranchAddress("Prompt_posX",&prompt_posx); 
        PromptTc->SetBranchAddress("Prompt_posY",&prompt_posy);
        PromptTc->SetBranchAddress("Prompt_posZAfterAVoffset",&prompt_posz); 
        PromptTc->SetBranchAddress("Prompt_eventid",&prompt_gtid);
        
        
        DelayTc->SetBranchAddress("Delay_energy",&delay_energy); 
        DelayTc->SetBranchAddress("delayedEcorr",&delayedEcorr);
        DelayTc->SetBranchAddress("Delay_clockCount50",&delay_clock50); 
        DelayTc->SetBranchAddress("Delay_posX",&delay_posx); 
        DelayTc->SetBranchAddress("Delay_posY",&delay_posy);
        DelayTc->SetBranchAddress("Delay_posZAfterAVoffset",&delay_posz);
        DelayTc->SetBranchAddress("Delay_eventid",&delay_gtid);
        DelayTc->SetBranchAddress("Delay_runID",&Delay_runID);
        
        
        int promptdatasize = PromptTc->GetEntries(); int delaydatasize = DelayTc->GetEntries();
        std::cout<<"Num of Tagged Events:  "<<PromptTc->GetEntries()<<std::endl;
        
        // for loop over lens of random number to select samples + fill into it
        
        int iEntry = (int)round(randomarr[i+1]*promptdatasize);
        int jEntry = (int)round(randomarr[i+2]*delaydatasize);
        PromptTc->GetEntry(iEntry); DelayTc->GetEntry(jEntry);
        // exclude the scenario of prompt and delay being same event
        if(prompt_gtid == delay_gtid) continue;
    
        //calculate dt
        dt = ((delay_clock50 - prompt_clock50) & 0x7FFFFFFFFFF)* 20.0;//ns
        promptR =  pow(( prompt_posx*prompt_posx + prompt_posy*prompt_posy + prompt_posz*prompt_posz),0.5 );
        //std::cout<<"prompt_posx "<<prompt_posx<<"prompt_posy "<<prompt_posy<<"prompt_posz "<<prompt_posz<<std::endl;
        delayR =  pow(( delay_posx*delay_posx + delay_posy*delay_posy + delay_posz*delay_posz),0.5 );
        dR = pow( (prompt_posx- delay_posx)*(prompt_posx- delay_posx)+ (prompt_posy- delay_posy)*(prompt_posy- delay_posy)+ (prompt_posz- delay_posz)*(prompt_posz- delay_posz) ,0.5);
        if( dR < 0.1) std::cout<<"prompt_gtid "<<prompt_gtid<<"delay_gtid "<<delay_gtid<<std::endl;

        dt_model = randomarr[i+3]*2E06;
        AccT->Fill();
    
        //PromptTc->Reset(); DelayTc->Reset();
        //std::cout<<"Num of Tagged Events:  "<<PromptTc->GetEntries()<<std::endl;
        inputFile->Close();
        //delete DelayTc;
        //delete PromptTc;
        
        
        
    }
    f_out->Write();
    f_out->Close();


}




int main(int argc, char** argv) {
    std::string inputrunlist = argv[1];
    std::string output_root_address = argv[2];
    std::string ratversion = argv[3];  // Old legacy entry, left for backwards compatibility with wrapper code (unused here)
    int gen_num = std::stoi(argv[4]); 
    
    
    // Addresses of simulation output files to be analysed
    std::vector<std::string> input_files;
    
    generate_accidental(inputrunlist,output_root_address,ratversion,gen_num);


    return 0;
}




