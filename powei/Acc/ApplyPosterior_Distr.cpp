
// generate posterior probability distribution for IBD and Accidental respectively
// apply posterier prob. to tree we loaded
// posterior prob = ln(L_reac/L_acc)+ ln(r_reac/r_acc) -> r_... is rate of ... at that run
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

float GetLLH(TH1D* IBDDelay_energyh1,TH2D* Analytical_IBDdRdth2,double delayedEcorr,double dR,double dt){
	Int_t bin_2D = Analytical_IBDdRdth2->FindBin(dR,dt,0);
	Int_t bin 	 = IBDDelay_energyh1->FindBin(delayedEcorr,0);
	//vstd::cout<<"bin "<<bin<<std::endl;
	std::cout<<dR<<" "<<dt<<std::endl;
	float Llh_dRdt   = Analytical_IBDdRdth2->GetBinContent(bin_2D);
	float Llh_energy = IBDDelay_energyh1->GetBinContent(bin);

	//std::cout<<"bin_2D "<<bin_2D <<std::endl;
	std::cout<<"Llh_dRdt"<<Llh_dRdt <<std::endl;
	std::cout<<"Llh_energy"<<Llh_energy <<std::endl;
	std::cout<<"log(Llh_dRdt*Llh_energy)"<<log(Llh_dRdt)+log(Llh_energy) <<std::endl;
	return log(Llh_dRdt*Llh_energy) ;
}
float GetAccLLH(TH1D* Acc_Delayenergyh1,TH2D* Acc_IBDdRdth2,double delayedEcorr,double dR,double dt_model){
	Int_t Accbin_2D = Acc_IBDdRdth2->FindBin(dR,dt_model,0);
	Int_t Accbin 	=  Acc_Delayenergyh1->FindBin(delayedEcorr,0);
	float AccLlh_dRdt   = Acc_IBDdRdth2->GetBinContent(Accbin_2D);
	float AccLlh_energy = Acc_Delayenergyh1->GetBinContent(Accbin);
	//std::cout<<"dR "<<dR <<std::endl;
	std::cout<<"AccLlh_dRdt"<<AccLlh_dRdt <<std::endl;
	std::cout<<"AccLlh_energy "<<AccLlh_energy<<std::endl;
	std::cout<<"log(AccLlh_dRdt*AccLlh_energy)"<<log(AccLlh_dRdt)+log(AccLlh_energy) <<std::endl;
	return log(AccLlh_dRdt*AccLlh_energy) ;
}
void FillBlank2D(TH2F* hist){
	float smallvalue = 0.1;
	for (int binX = 1; binX <= hist->GetNbinsX(); ++binX) {
        for (int binY = 1; binY <= hist->GetNbinsY(); ++binY) {
			//std::cout<<binX<<" "<<binY<<std::endl;
            hist->AddBinContent(hist->GetBin(binX, binY), smallvalue);
        }
    }

}
void GetContent2D(TH2F* hist){
	float smallvalue = 0.1;
	for (int binX = 1; binX <= hist->GetNbinsX(); ++binX) {
        for (int binY = 1; binY <= hist->GetNbinsY(); ++binY) {
			std::cout<<"GetBinContent "<<hist->GetBinContent(hist->GetBin(binX, binY))<<std::endl;
        }
    }

}

void FillBlank1D(TH1F* hist){
	float smallvalue = 0.1;
	for (int binX = 1; binX <= hist->GetNbinsX(); ++binX) {
        hist->AddBinContent(hist->GetBin(binX), smallvalue);
    }
    

}

std::vector<std::pair<int, double>> GetRunID_AccIbdRatio(std::string listdir){
    
	std::vector<std::pair<int, double>> RunID_AccIbdRatio;
	std::ifstream infile(listdir);

    if (!infile) {
        std::cerr << "Unable to open file: " << listdir << std::endl;
        exit(0); // Exit with error code
    }

    std::string line;
    while (std::getline(infile, line)) {
        std::istringstream iss(line);
        std::string Runid; std::string AccIbdRatio;

        while (iss >> Runid >>AccIbdRatio) {
			std::pair<int, double> Pair(std::stoi(Runid), std::stod(AccIbdRatio));
            RunID_AccIbdRatio.push_back(Pair);
        }

        // Process or print the elements vector
		/*
        std::cout << "Line elements: ";
        for (const auto& el : RunID_AccIbdRatio) {
            std::cout << el.first << " ";
        }
        std::cout << std::endl;
		*/
    }

    infile.close();
    return RunID_AccIbdRatio; // Successful execution
}
double findRateRatiofromRunId(std::vector<std::pair<int, double>> RunID_AccIbdRatio,int RunID){
    
    int target = RunID;
    auto it = std::find_if(RunID_AccIbdRatio.begin(), RunID_AccIbdRatio.end(),
        [target](const std::pair<int, double>& pairelments) {
            return pairelments.first == target;
                           });
    if (it != RunID_AccIbdRatio.end()) {
        return it->second; // Return the corresponding double value
    } else {
        std::cout<< "not found ratio in Run "<< RunID<<std::endl; // Return std::nullopt if the element is not found
		return 0.;
	}
}

void ClassifierIBD(std::string inputfilepath, std::string output_root_address, std::string Proc_rat){
	// set the output file, tree
	TFile *f_out = new TFile( output_root_address.c_str(),"RECREATE");
	// Prompt_Tc and Delay_Tc are Tree of inputfile argument(can be sim or data)
	TChain *Prompt_Tc   = new TChain("PromptT");
    TChain *Delay_Tc    = new TChain("DelayT");
	Prompt_Tc->Add((inputfilepath).c_str());
	Delay_Tc->Add((inputfilepath).c_str());
	TTree *NewPrompt_Tc = Prompt_Tc->CloneTree(0);
	TTree *NewDelay_Tc  = Delay_Tc->CloneTree(0);
	// NewPrompt_Tc and NewDelay_Tc are output Tree
	TFile *oldfile = new TFile(inputfilepath.c_str());
    TTree *NewPrompt_Tc_test = (TTree*)oldfile->Get("PromptT")->Clone(0); 
	TTree *NewDelay_Tc_test = (TTree*)oldfile->Get("DelayT")->Clone(0);
	if (Prompt_Tc->GetEntries()== 0 || Delay_Tc->GetEntries()== 0){
		std::cout<<"No Entries in this file... no need to calulate posterior prob..."<<std::endl;
   		f_out->cd();
		//NewPrompt_Tc->SetName("PromptT");
		//NewDelay_Tc->SetName("DelayT");
		NewPrompt_Tc_test->Write();
		NewDelay_Tc_test->Write();
   		//f_out->Write();
    	f_out->Close();
		return;
	}

	double dR; double dt; int RunID; double delayedEcorr;
	
	Delay_Tc->SetBranchAddress("dt",&dt); 
	Delay_Tc->SetBranchAddress("dR",&dR); ; 
	Delay_Tc->SetBranchAddress("delayedEcorr",&delayedEcorr); 
	Delay_Tc->SetBranchAddress("Delay_runID",&RunID);
	
	std::string IBDinputpath = "/data/snoplus3/weiii/antinu/mycuts/Ntuple_sim/Reactoribd";
	std::string Accinputpath = "/data/snoplus3/weiii/antinu/mycuts/Ntuple_Acc/accidental";
	
	
	float Log_Llh_IBD; float Log_Llh_Acc; float post_prob;
	double posterior; double likelihood_Acc; double likelihood_IBD;

    TFile *file_delayE = TFile::Open("/home/huangp/AntiNu/Analysis/delayedPDFs.root");
	TH1D  *Acc_Delayenergyh1      = new TH1D();
	TH1D *IBDDelay_energyh1  = new TH1D();
	// load histo
	file_delayE->GetObject("Acc_Delayenergyh1", Acc_Delayenergyh1);
	file_delayE->GetObject("IBDDelay_energyh1", IBDDelay_energyh1);
	//std::cout<<"Entries of Acc_Delayenergyh1 "<<Acc_Delayenergyh1->GetEntries()<<std::endl;

	TFile *file_dRdt = TFile::Open("/home/huangp/AntiNu/Analysis/dtdrpdf.root");
	TH2D *Analytical_IBDdRdth2 = new TH2D();
	TH2D *Acc_IBDdRdth2        = new TH2D();
	// load histo
	file_dRdt->GetObject("Analytical_IBDdRdth2", Analytical_IBDdRdth2);
	file_dRdt->GetObject("Acc_IBDdRdth2", Acc_IBDdRdth2);



	


	// Branch for saving
	//NewDelay_Tc->Branch("dR",&dR);
    //NewDelay_Tc->Branch("dt_model",&dt_model);
    //NewDelay_Tc->Branch("delayedEcorr",&delayE);
    //NewDelay_Tc->Branch("RunID",&RunID);
	NewDelay_Tc->Branch("Log_Llh_IBD",&Log_Llh_IBD);
	NewDelay_Tc->Branch("Log_Llh_Acc",&Log_Llh_Acc);
	NewDelay_Tc->Branch("post_prob",&post_prob);


	//NewPrompt_Tc->Branch("dR",&dR);
    //NewPrompt_Tc->Branch("dt",&dt);
    //NewPrompt_Tc->Branch("delayE",&delayE);
    //NewPrompt_Tc->Branch("RunID",&RunID);
	NewPrompt_Tc->Branch("Log_Llh_IBD",&Log_Llh_IBD);
	NewPrompt_Tc->Branch("Log_Llh_Acc",&Log_Llh_Acc);
	NewPrompt_Tc->Branch("post_prob",&post_prob);

	
	// load reactorIBD/acc rate
	std::vector<std::pair<int, double>> RunID_AccIbdRatio;
	RunID_AccIbdRatio = GetRunID_AccIbdRatio("/home/huangp/AntiNu/Analysis/RateRatio.txt");

	// load rate of accidental and IBD rate with RUN_NUM
	std::cout<< "Ready to Processing IBD Samples"<<std::endl;
	for(int iEntry = 0; iEntry < Delay_Tc->GetEntries(); iEntry++){
		//std::cout<<iEntry<<std::endl;
		
		Delay_Tc->GetEntry(iEntry);
		post_prob = 0;
		std::cout<<delayedEcorr<<std::endl;
		if(delayedEcorr>2.5 || delayedEcorr<1.85) post_prob = -999999;
		if(dR>2500) post_prob = -999999;
		if(dt>2000000) post_prob = -999999;
		if(post_prob != -999999){
			Log_Llh_IBD = GetLLH(IBDDelay_energyh1,Analytical_IBDdRdth2,delayedEcorr,dR,dt);
			Log_Llh_Acc = GetAccLLH(Acc_Delayenergyh1,Acc_IBDdRdth2,delayedEcorr,dR,dt);
			//double IBD_ACCRateRatio = findRateRatiofromRunId(RunID_AccIbdRatio,RunID);
			//if (IBD_ACCRateRatio - 0. < 0.0001){
				//post_prob = Log_Llh_IBD - Log_Llh_Acc - log(IBD_ACCRateRatio);
			//}
			//else{
			post_prob = Log_Llh_IBD - Log_Llh_Acc;
			//}
		}
		std::cout<<post_prob<<std::endl;
		std::cout<<"Progress of IBD "<<float(iEntry)/float(Delay_Tc->GetEntries())<<std::endl;
		NewDelay_Tc->Fill();
		Prompt_Tc->GetEntry(iEntry);
		NewPrompt_Tc->Fill();
		
	}
	f_out->cd();
	NewPrompt_Tc->SetName("PromptT");
	NewDelay_Tc->SetName("DelayT");
	NewDelay_Tc->Write();
	NewPrompt_Tc->Write();
	//f_out->Write();
    f_out->Close();

	delete f_out;
	//delete Prompt_Tc;
	//delete Delay_Tc;
	//delete NewPrompt_Tc;
	//delete NewDelay_Tc;
	delete Acc_Delayenergyh1;
	delete Acc_IBDdRdth2;
	delete IBDDelay_energyh1;
	delete Analytical_IBDdRdth2;
	
 
	
}

void ClassifierAcc(std::string inputfilepath, std::string output_root_address, std::string Proc_rat){
	// set the output file, tree
	TFile *f_out = new TFile( output_root_address.c_str(),"RECREATE");
	
	// Acc_Tc is Tree of inputfile argument(only applied to accidental)
    TChain *Acc_Tc    = new TChain("AccT");
	Acc_Tc->Add((inputfilepath).c_str());
	TTree *NewAcc_Tc  = Acc_Tc->CloneTree(0);
	// NewAcc_Tc is output Tree
	TFile *oldfile = new TFile(inputfilepath.c_str());
	TTree *NewAcc_Tc_test = (TTree*)oldfile->Get("AccT")->Clone(0);
	
	if (Acc_Tc->GetEntries()== 0){
		std::cout<<"No Entries in this file... no need to calulate posterior prob..."<<std::endl;
   		f_out->cd();
		//NewDelay_Tc->SetName("AccT");
		NewAcc_Tc_test->Write();
   		//f_out->Write();
    	f_out->Close();
		return;
	}
	
	
	double dR; double dt_model; int RunID; double delayedEcorr;
	
	Acc_Tc->SetBranchAddress("dt_model",&dt_model); 
	Acc_Tc->SetBranchAddress("dR",&dR); ; 
	Acc_Tc->SetBranchAddress("delayedEcorr",&delayedEcorr); 
	Acc_Tc->SetBranchAddress("RunID",&RunID);
	
	std::string IBDinputpath = "/data/snoplus3/weiii/antinu/mycuts/Ntuple_sim/Reactoribd";
	std::string Accinputpath = "/data/snoplus3/weiii/antinu/mycuts/Ntuple_Acc/accidental";
	
	
	float Log_Llh_IBD; float Log_Llh_Acc; float post_prob;
	double posterior; double likelihood_Acc; double likelihood_IBD;

    TFile *file_delayE = TFile::Open("/home/huangp/AntiNu/Analysis/delayedPDFs.root");
	TH1D  *Acc_Delayenergyh1      = new TH1D();
	TH1D *IBDDelay_energyh1  = new TH1D();
	// load histo
	file_delayE->GetObject("Acc_Delayenergyh1", Acc_Delayenergyh1);
	file_delayE->GetObject("IBDDelay_energyh1", IBDDelay_energyh1);
	//std::cout<<"Entries of Acc_Delayenergyh1 "<<Acc_Delayenergyh1->GetEntries()<<std::endl;

	TFile *file_dRdt = TFile::Open("/home/huangp/AntiNu/Analysis/dtdrpdf.root");
	TH2D *Analytical_IBDdRdth2 = new TH2D();
	TH2D *Acc_IBDdRdth2        = new TH2D();
	// load histo
	file_dRdt->GetObject("Analytical_IBDdRdth2", Analytical_IBDdRdth2);
	file_dRdt->GetObject("Acc_IBDdRdth2", Acc_IBDdRdth2);



	


	// Branch for saving
	//NewAcc_Tc->Branch("dR",&dR);
    //NewAcc_Tc->Branch("dt_model",&dt_model);
    //NewAcc_Tc->Branch("delayedEcorr",&delayE);
    //NewAcc_Tc->Branch("RunID",&RunID);
	NewAcc_Tc->Branch("Log_Llh_IBD",&Log_Llh_IBD);
	NewAcc_Tc->Branch("Log_Llh_Acc",&Log_Llh_Acc);
	NewAcc_Tc->Branch("post_prob",&post_prob);


	

	
	// load reactorIBD/acc rate
	std::vector<std::pair<int, double>> RunID_AccIbdRatio;
	RunID_AccIbdRatio = GetRunID_AccIbdRatio("/home/huangp/AntiNu/Analysis/RateRatio.txt");
	
	// load rate of accidental and IBD rate with RUN_NUM
	std::cout<< "Ready to Processing IBD Samples"<<std::endl;
	for(int iEntry = 0; iEntry < Acc_Tc->GetEntries(); iEntry++){
		//std::cout<<iEntry<<std::endl;
		post_prob = 0;
		Acc_Tc->GetEntry(iEntry);
		std::cout<<delayedEcorr<<std::endl;
		if(delayedEcorr>2.5 || delayedEcorr<1.85) post_prob = -999999;
		if(dR>2500) post_prob = -999999;
		if(dt_model>2000000) post_prob = -999999;
		if(post_prob != -999999){
			Log_Llh_IBD = GetLLH(IBDDelay_energyh1,Analytical_IBDdRdth2,delayedEcorr,dR,dt_model);
			Log_Llh_Acc = GetAccLLH(Acc_Delayenergyh1,Acc_IBDdRdth2,delayedEcorr,dR,dt_model);
			//double IBD_ACCRateRatio = findRateRatiofromRunId(RunID_AccIbdRatio,RunID);
			//if (IBD_ACCRateRatio - 0. < 0.0001){
				//post_prob = Log_Llh_IBD - Log_Llh_Acc - log(IBD_ACCRateRatio);
			//}
			//else{
			post_prob = Log_Llh_IBD - Log_Llh_Acc;
			//}
		}
		std::cout<<post_prob<<std::endl;
		std::cout<<"Progress of IBD "<<float(iEntry)/float(Acc_Tc->GetEntries())<<std::endl;
		NewAcc_Tc->Fill();
		
	}
	
	f_out->cd();
	NewAcc_Tc->SetName("AccT");
	NewAcc_Tc->Write();
	//f_out->Write();
    f_out->Close();
	
	//delete f_out;
	//delete Acc_Tc;
	//delete NewAcc_Tc;
	
	delete Acc_Delayenergyh1;
	delete Acc_IBDdRdth2;
	delete IBDDelay_energyh1;
	delete Analytical_IBDdRdth2;
 
	
}



int main(int argc, char** argv) {
	std::string inputfilepath = argv[1];
	std::string output_root_address = argv[2];
	std::string ratversion = argv[3];  // Old legacy entry, left for backwards compatibility with wrapper code (unused here)
	std::string eventtype  = argv[4]; //if eventtype == "Acc" , using different branch name to get entry inside
	
	// Addresses of simulation output files to be analysed
	std::vector<std::string> input_files;
	std::cout<<"eventtype "<<eventtype<<std::endl;
	if (eventtype == "Acc"){
		ClassifierAcc(inputfilepath,output_root_address,ratversion);
	}
	else {
		ClassifierIBD(inputfilepath,output_root_address,ratversion);
	}


	return 0;
}
