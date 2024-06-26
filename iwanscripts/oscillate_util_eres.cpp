#include <TFile.h>
#include <TMath.h>
#include <string>
#include <TNtuple.h>
#include <TTree.h>
#include <iostream>
#include <TObject.h>
//#include <CLHEP/Random/Randomize.h>
#include <TRandom3.h>

Double_t correction_to_EPrompt (Double_t x){
  //double correc = 8.83411e-02 + (9.02287e-03)*x;   //vs EPrompt
  double correc = 0.;
  return correc;
}

Double_t correction_to_truth (Double_t x){
  //double correc = (-1.)*(7.63260e-02 + (1.19968e-02)*x);   //vs antinuKE - 0.784
  double correc = (-1.)*(0.0638836 + (-0.0200583)*x);   //vs positronKE - 1.022 //v6.17.6
  //double correc = 0.;
  return correc;
}

Double_t NuSurvProb(Double_t nuE, Double_t baseline, Double_t del_m_sqr_21, Double_t sin_sqr_theta_12, Double_t sin_sqr_theta_13){
    Double_t f_s_sqr2_theta_12 = pow(sin(2.0 * TMath::ASin(sqrt(sin_sqr_theta_12))), 2.0);
    Double_t f_s4 = pow(sin_sqr_theta_13, 2.0);
    Double_t f_c4 = pow(1.0-sin_sqr_theta_13, 2.0);
    Double_t scale = 1.267e3; // for nuE in [MeV] and baseline in [km]
    Double_t s_sqr_dm_be = pow(sin(scale * del_m_sqr_21 * baseline / nuE), 2.0);
    Double_t f_osc_prob = (f_c4 * (1.0 - f_s_sqr2_theta_12 * s_sqr_dm_be) + f_s4);
    return f_osc_prob;
}

Double_t NuSurvProb_geo(Double_t sin_sqr_theta_12){
    Double_t f_s_sqr2_theta_12 = pow(sin(2.0 * TMath::ASin(sqrt(sin_sqr_theta_12))), 2.0);
    Double_t f_osc_prob = 1.0 - (0.5*f_s_sqr2_theta_12);
    return f_osc_prob;
}

void ntOscillate_pruned(TTree *in_tree, TNtuple *out_tree_prompt, Double_t del_m_sqr_21, Double_t sin_sqr_theta_12, Double_t sin_sqr_theta_13, Double_t fixed_distance, const bool apply_energy_resolution_convolution, const Double_t annihilation_energy) {

    //
    // takes a TTree with multiple branches, oscillates using KE branch, fills two TNtuples each with a single branch
    //
    ULong64_t n_entries = in_tree->GetEntries();
    Double_t surv_prob, mc_energy_nu, mc_energy_positron, ev_energy_p1, distance;
    in_tree->SetBranchAddress("mc_neutrino_energy", &mc_energy_nu);
    in_tree->SetBranchAddress("mc_positron_energy", &mc_energy_positron);
    in_tree->SetBranchAddress("ev_fit_energy_p1", &ev_energy_p1);

    if (fixed_distance>0)
        distance = fixed_distance;
    else
        in_tree->SetBranchAddress("reactor_info_distance", &distance);

    TRandom3 *random_generator = new TRandom3();
    random_generator->SetSeed(0);
    
    for (ULong64_t i = 0; i < n_entries; i++){
        in_tree->GetEntry(i);
        surv_prob = NuSurvProb(mc_energy_nu, distance, del_m_sqr_21, sin_sqr_theta_12, sin_sqr_theta_13);

        //const Double_t random = CLHEP::HepUniformRand();
        Double_t random = random_generator->Uniform();
	
        if (surv_prob > random){
            //printf("ke:%0.5f ev:%0.5f diff:%0.5f\n", (Float_t)mc_energy_nu, (Float_t)ev_energy_p1, (Float_t)mc_energy_nu-(Float_t)ev_energy_p1);
            //out_tree_ke->Fill((Float_t)mc_energy_nu);
	    if (apply_energy_resolution_convolution){
	        //Double_t mc_energy_p1 = mc_energy_nu + annihilation_energy;
	        Double_t mc_energy_p1 = mc_energy_positron + annihilation_energy;
		//std::cout<<"oscillate: "<<(mc_energy_p1 + correction_to_truth(mc_energy_p1))<<" <->"<<ev_energy_p1<<std::endl;
out_tree_prompt->Fill((Float_t)(mc_energy_p1 + correction_to_truth(mc_energy_positron)));
	    }else{
	        out_tree_prompt->Fill((Float_t)ev_energy_p1);
	    }
        }
    }
}

void ntOscillate_pruned(TTree *in_tree, TNtuple *out_tree_prompt, Double_t del_m_sqr_21, Double_t sin_sqr_theta_12, Double_t sin_sqr_theta_13, Double_t fixed_distance) {
    ntOscillate_pruned(in_tree, out_tree_prompt, del_m_sqr_21, sin_sqr_theta_12, sin_sqr_theta_13, fixed_distance, false, 1.022);
}

void ntOscillate_pruned(TTree *in_tree, TNtuple *out_tree_prompt, Double_t del_m_sqr_21, Double_t sin_sqr_theta_12, Double_t sin_sqr_theta_13) {
    ntOscillate_pruned(in_tree, out_tree_prompt, del_m_sqr_21, sin_sqr_theta_12, sin_sqr_theta_13, -9000, false, 1.022);
}

void ntOscillate_pruned_geo(TTree *in_tree, TNtuple *out_tree_prompt, Double_t sin_sqr_theta_12, const bool apply_energy_resolution_convolution, const Double_t annihilation_energy) {

    ULong64_t n_entries = in_tree->GetEntries();
    Double_t surv_prob, mc_energy_nu, mc_energy_positron, ev_energy_p1;
    in_tree->SetBranchAddress("mc_neutrino_energy", &mc_energy_nu);
    in_tree->SetBranchAddress("mc_positron_energy", &mc_energy_positron);
    in_tree->SetBranchAddress("ev_fit_energy_p1", &ev_energy_p1);
    
    TRandom3 *random_generator = new TRandom3();
    random_generator->SetSeed(0);

    for (ULong64_t i = 0; i < n_entries; i++){
        in_tree->GetEntry(i);

	surv_prob = NuSurvProb_geo(sin_sqr_theta_12);  
        Double_t random = random_generator->Uniform();

        if (surv_prob > random){
	    if (apply_energy_resolution_convolution){
	        Double_t mc_energy_p1 = mc_energy_positron + annihilation_energy;
	        out_tree_prompt->Fill((Float_t)(mc_energy_p1 + correction_to_truth(mc_energy_positron)));
	    }else{
	        out_tree_prompt->Fill((Float_t)ev_energy_p1);
	    }
	}
    }
}

void ntOscillate_pruned_geo(TTree *in_tree, TNtuple *out_tree_prompt, Double_t sin_sqr_theta_12) {
    ntOscillate_pruned_geo(in_tree, out_tree_prompt, sin_sqr_theta_12, false, 1.022);
}

void ntNoOscillate_pruned(TTree *in_tree, TNtuple *out_tree_prompt, const bool apply_energy_resolution_convolution, const Double_t annihilation_energy) {
    //
    // takes a TTree with multiple branches, doesn't apply any oscillation, to made osc_pdf == unosc_pdf, keep things tidy in MultiOscillation.cpp
    //
    ULong64_t n_entries = in_tree->GetEntries();
    Double_t surv_prob, mc_energy_nu, mc_energy_positron, ev_energy_p1;
    in_tree->SetBranchAddress("mc_neutrino_energy", &mc_energy_nu);
    in_tree->SetBranchAddress("mc_positron_energy", &mc_energy_positron);
    in_tree->SetBranchAddress("ev_fit_energy_p1", &ev_energy_p1);
    
    for (ULong64_t i = 0; i < n_entries; i++){
        in_tree->GetEntry(i);
	if (apply_energy_resolution_convolution){
	    Double_t mc_energy_p1 = mc_energy_nu + annihilation_energy;
	    out_tree_prompt->Fill((Float_t)(mc_energy_p1 + correction_to_truth(mc_energy_positron)));
	}else{
	    out_tree_prompt->Fill((Float_t)ev_energy_p1);
	}
    }
}

void ntNoOscillate_pruned(TTree *in_tree, TNtuple *out_tree_prompt) {
    ntNoOscillate_pruned(in_tree, out_tree_prompt, false, 1.022);
}

void ntOscillate(TTree *in_tree, TTree *out_tree, Double_t del_m_sqr_21, Double_t sin_sqr_theta_12, Double_t sin_sqr_theta_13, Double_t fixed_distance) {

    ULong64_t n_entries = in_tree->GetEntries();
    Double_t surv_prob, mc_energy_nu, distance;
    in_tree->SetBranchAddress("mc_neutrino_energy", &mc_energy_nu);

    if (fixed_distance>0)
        distance = fixed_distance;
    else
        in_tree->SetBranchAddress("reactor_info_distance", &distance);
    
    TRandom3 *random_generator = new TRandom3();
    random_generator->SetSeed(0);
    
    for (ULong64_t i = 0; i < n_entries; i++){
        in_tree->GetEntry(i);

	surv_prob = NuSurvProb(mc_energy_nu, distance, del_m_sqr_21, sin_sqr_theta_12, sin_sqr_theta_13);  
        //const Double_t random = CLHEP::HepUniformRand();
        Double_t random = random_generator->Uniform();

        if (surv_prob > random)
            out_tree->Fill();
    }
}

void ntOscillate_pdf(TTree *in_tree, TTree *out_tree, Double_t del_m_sqr_21, Double_t sin_sqr_theta_12, Double_t sin_sqr_theta_13) {

    ULong64_t n_entries = in_tree->GetEntries();
    Double_t surv_prob, mc_energy_nu, distance;
    in_tree->SetBranchAddress("mc_neutrino_energy", &mc_energy_nu);

    in_tree->SetBranchAddress("reactor_info_distance", &distance);
    
    TRandom3 *random_generator = new TRandom3();
    random_generator->SetSeed(0);
    
    for (ULong64_t i = 0; i < n_entries; i++){
        in_tree->GetEntry(i);

	if (distance > 0)
	  surv_prob = NuSurvProb(mc_energy_nu, distance, del_m_sqr_21, sin_sqr_theta_12, sin_sqr_theta_13);
	else
	  surv_prob = NuSurvProb_geo(sin_sqr_theta_12);
	  
        //const Double_t random = CLHEP::HepUniformRand();
        Double_t random = random_generator->Uniform();

        if (surv_prob > random)
            out_tree->Fill();
    }
}

void ntOscillate_geo(TTree *in_tree, TTree *out_tree, Double_t sin_sqr_theta_12) {

    ULong64_t n_entries = in_tree->GetEntries();

    TRandom3 *random_generator = new TRandom3();
    random_generator->SetSeed(0);

    Double_t surv_prob = NuSurvProb_geo(sin_sqr_theta_12);

    for (ULong64_t i = 0; i < n_entries; i++){
        in_tree->GetEntry(i);
	Double_t random = random_generator->Uniform();

        if (surv_prob > random)
	    out_tree->Fill();
    }
}

void ntOscillate(TTree *in_tree, TTree *out_tree, Double_t del_m_sqr_21, Double_t sin_sqr_theta_12, Double_t sin_sqr_theta_13) {
    ntOscillate(in_tree, out_tree, del_m_sqr_21, sin_sqr_theta_12, sin_sqr_theta_13, -9000);
}

void write_file(const char* nt_in, const char* nt_out, Double_t del_m_sqr_21, Double_t sin_sqr_theta_12, Double_t sin_sqr_theta_13) {

    TFile *f_in = new TFile(nt_in);
    TTree *in_tree = (TTree*)f_in->Get("nt");
    TFile *f_out = new TFile(nt_out, "RECREATE");
    TTree *out_tree = in_tree->CloneTree(0);

    ntOscillate(in_tree, out_tree, del_m_sqr_21, sin_sqr_theta_12, sin_sqr_theta_13);

    f_out->cd();
    out_tree->Write();
    f_in->Close();
    f_out->Close();
    delete f_in;
    delete f_out;
}

void write_file(const char* nt_in, const char* nt_out, Double_t del_m_sqr_21, Double_t sin_sqr_theta_12, Double_t sin_sqr_theta_13, Double_t fixed_distance) {

    TFile *f_in = new TFile(nt_in);
    TTree *in_tree = (TTree*)f_in->Get("nt");
    TFile *f_out = new TFile(nt_out, "RECREATE");
    TTree *out_tree = in_tree->CloneTree(0);

    ntOscillate(in_tree, out_tree, del_m_sqr_21, sin_sqr_theta_12, sin_sqr_theta_13, fixed_distance);

    f_out->cd();
    out_tree->Write();
    f_in->Close();
    f_out->Close();
    delete f_in;
    delete f_out;
}

void write_file_pdf(const char* nt_in, const char* nt_out, Double_t del_m_sqr_21, Double_t sin_sqr_theta_12, Double_t sin_sqr_theta_13) {

    TFile *f_in = new TFile(nt_in);
    TTree *in_tree = (TTree*)f_in->Get("nt");
    TFile *f_out = new TFile(nt_out, "RECREATE");
    TTree *out_tree = in_tree->CloneTree(0);

    ntOscillate_pdf(in_tree, out_tree, del_m_sqr_21, sin_sqr_theta_12, sin_sqr_theta_13);

    f_out->cd();
    out_tree->Write();
    f_in->Close();
    f_out->Close();
    delete f_in;
    delete f_out;
}

void write_file_geo(const char* nt_in, const char* nt_out, Double_t sin_sqr_theta_12) {

    TFile *f_in = new TFile(nt_in);
    TTree *in_tree = (TTree*)f_in->Get("nt");
    TFile *f_out = new TFile(nt_out, "RECREATE");
    TTree *out_tree = in_tree->CloneTree(0);

    ntOscillate_geo(in_tree, out_tree, sin_sqr_theta_12);

    f_out->cd();
    out_tree->Write();
    f_in->Close();
    f_out->Close();
    delete f_in;
    delete f_out;
}

void write_file_pruned(const char* nt_in, const char* nt_prompt_out, Double_t del_m_sqr_21, Double_t sin_sqr_theta_12, Double_t sin_sqr_theta_13, Double_t distance) {

    TFile *f_in = new TFile(nt_in);
    TTree *in_tree = (TTree*)f_in->Get("nt");

    // the ordering of things here to keep ROOT happy is exactly the following (TFile then TTree, TFile again TTree)
    //TFile *f_ke_out = new TFile(nt_ke_out, "RECREATE");
    //TNtuple *out_tree_ke = new TNtuple("nt","Anti-neutrino processed tree", "mc_neutrino_energy");

    TFile *f_prompt_out = new TFile(nt_prompt_out, "RECREATE");
    TNtuple *out_tree_prompt = new TNtuple("nt", "Oscillated Prompt Energy", "ev_fit_energy_p1");

    if (distance<0)
        ntOscillate_pruned(in_tree, out_tree_prompt, del_m_sqr_21, sin_sqr_theta_12, sin_sqr_theta_13);
    else
        ntOscillate_pruned(in_tree, out_tree_prompt, del_m_sqr_21, sin_sqr_theta_12, sin_sqr_theta_13, distance, false, 1.022);

    //f_ke_out->cd();
    //out_tree_ke->Write();
    //f_ke_out->Close();
    //delete f_ke_out;

    f_prompt_out->cd();
    out_tree_prompt->Write();
    f_prompt_out->Close();
    delete f_prompt_out;

    f_in->cd();
    f_in->Close();
    delete f_in;
}

void write_file_pruned(const char* nt_in, const char* nt_prompt_out, Double_t del_m_sqr_21, Double_t sin_sqr_theta_12, Double_t sin_sqr_theta_13) {
    write_file_pruned(nt_in, nt_prompt_out, del_m_sqr_21, sin_sqr_theta_12, sin_sqr_theta_13, -9000);
}
