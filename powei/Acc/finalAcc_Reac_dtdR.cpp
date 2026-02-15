#include <TROOT.h>
#include <TChain.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TF1.h>
#include <TFile.h>
#include <iostream>
#include <vector>
#include <glob.h>

void Fill(TH2F* h2, float offset){
	int binx = h2->GetNbinsX();  int biny = h2->GetNbinsY();
	//std::cout<<"h2binsX"<<binx<<" "<<biny<<std::endl;
	for(int x = 1; x<= binx; x++){
		for(int y = 1; y<= biny; y++){
			Int_t gbin = h2->GetBin(x,y);
			float val = h2->GetBinContent(gbin);
			h2->SetBinContent(gbin,val+offset);
		}
	}
}
// Function to load simulation files
void loadfile(TChain* treename, const std::string& fileloc) {
    glob_t glob_result;
    std::string pattern = fileloc + "scaled_oscillated*.root";
    glob(pattern.c_str(), GLOB_TILDE, NULL, &glob_result);
    
    for(unsigned int i = 0; i < glob_result.gl_pathc; ++i) {
        std::string f = glob_result.gl_pathv[i];
        if(std::stoi(f.substr(f.size() - 24, 6)) > 300600) {
            treename->Add(f.c_str());
        }
    }
    globfree(&glob_result);
}

// Function to load data files
void loaddatafile(TChain* treename, const std::string& fileloc) {
    glob_t glob_result;
    std::string pattern = fileloc + "*.root";
    glob(pattern.c_str(), GLOB_TILDE, NULL, &glob_result);

    for(unsigned int i = 0; i < glob_result.gl_pathc; ++i) {
        std::string f = glob_result.gl_pathv[i];
        treename->Add(f.c_str());
    }
    globfree(&glob_result);
}


Double_t expdecay_func(Double_t t, Double_t par)
{
   Float_t tt =t;
   Double_t f = par*exp(-tt/214150.);
   return f;
}
void FillIBD_2dhistocolumn(TH2D* h2,TF1* fit_func,double A, int ibinX){
    //x:dR[mm], y:dt[ns]
    int biny = h2->GetNbinsY();
    
    for(int y = 1; y<= biny; y++){
        Int_t gbin = h2->GetBin(ibinX,y);
        //float val = h2->GetBinContent(gbin);
        Double_t bcx = ((TAxis*)h2->GetXaxis())->GetBinCenter(ibinX);
        Double_t bcy = ((TAxis*)h2->GetYaxis())->GetBinCenter(y);
        Double_t llh = expdecay_func(bcy, A);
        //std::cout<<bcx<<" "<<bcy<<std::endl;
        //std::cout<<"llh "<<llh<<std::endl;
        h2->SetBinContent(gbin,llh);
    }
	
}
void FillAcc_2dhistocolumn(TH2D* h2,TF1* fit_func,double A, int ibinX){
    //x:dR[mm], y:dt[ns]
    int biny = h2->GetNbinsY();
    
    for(int y = 1; y<= biny; y++){
        Int_t gbin = h2->GetBin(ibinX,y);
        //float val = h2->GetBinContent(gbin);
        Double_t llh = A;
        //std::cout<<bcx<<" "<<bcy<<std::endl;
        //std::cout<<"llh "<<llh<<std::endl;
        h2->SetBinContent(gbin,llh);
    }
	
}
/*
void addoffset(TH2F* h2, float offset){
	int binx = h2->GetNbinsX();  int biny = h2->GetNbinsY();
	//std::cout<<"h2binsX"<<binx<<" "<<biny<<std::endl;
	for(int x = 1; x<= binx; x++){
		for(int y = 1; y<= biny; y++){
			Int_t gbin = h2->GetBin(x,y);
			float val = h2->GetBinContent(gbin);
            
			h2->SetBinContent(gbin,val+offset);
		}
	}
}
*/
void finalAcc_Reac_dtdR() {
    // Event types
    gStyle->SetOptStat(0);
    std::vector<std::string> antinueventtype = {
        "Geoibd_U", "Geoibd_Th", "Tl210", "Alphan_Lab_13c",
        "Alphan_Lab_Avout_Av_18o", "Alphan_Lab_Avin_Av_13c",
        "Alphan_Lab_Avout_Av_13c", "Alphan_Lab_Avin_Av_18o",
        "Reactoribd", "Bipo214", "Bipo212"
    };

    std::string filepath = "/data/snoplus3/weiii/antinu/mycuts/Ntuple_sim";
    std::string datafilepath = "/data/snoplus3/weiii/antinu/mycuts/Ntuple_Acc/accidental/";

    // TChain for Reactor IBD
    TChain IBDDelayT("PromptT");
    for (const auto& itype : antinueventtype) {
        std::string fileloc = filepath + "/" + itype + "/";
        if (itype == "Reactoribd") {
            loadfile(&IBDDelayT, fileloc);
        }
    }

    // TChain for Accidental Data
    TChain DataDelayT("AccT");
    loaddatafile(&DataDelayT, datafilepath);
    
    // Histograms for dt distribution
    int dRbin = 100, dtbin = 100;
    
    // Canvas 2: dR plot
    TCanvas c2("c2", "", 800, 500);
    TH1F IBDdRh1("IBDdRh1", "", dRbin, 0., 2500);
    TH1F Data_dRh1("Data_dRh1", "", dRbin, 0., 2500);

    IBDDelayT.Project("IBDdRh1", "dR");
    DataDelayT.Project("Data_dRh1", "dR", "");

    //IBDdRh1.Scale(1. / IBDdRh1.Integral(), "nosw2");
    //Data_dRh1.Scale(1. / Data_dRh1.Integral(), "nosw2");

    IBDdRh1.SetLineWidth(2);
    Data_dRh1.SetLineWidth(2);
    IBDdRh1.SetLineColor(1);
    Data_dRh1.SetLineColor(2);
    IBDdRh1.GetXaxis()->SetTitle("Delta R [mm]");
    IBDdRh1.GetYaxis()->SetTitle("Asimov Events/10 [mm]");

    IBDdRh1.Draw();
    Data_dRh1.Draw("same");

    TLegend legend1(0.65, 0.65, 0.85, 0.85);
    legend1.AddEntry(&IBDdRh1, "Reactor antinu", "L");
    legend1.AddEntry(&Data_dRh1, "Data", "L");
    legend1.Draw();

    //c2.SaveAs("/home/huangp/AntiNu/Plots/Acc_Reac_dR.png");

    // Canvas 3: 2D plot Reactor dt-dR
    TCanvas c3("c3", "", 800, 500);
    c3.SetLogz();
    TH2D IBDdRdth2("IBDdRdth2", "", dRbin, 0., 2500, dtbin, 0., 2E6);

    IBDDelayT.Project("IBDdRdth2", "dt:dR");
    
    IBDdRdth2.Scale(1. / IBDdRdth2.Integral(), "nosw2");
    IBDdRdth2.GetXaxis()->SetTitle("Delta R [mm]");
    IBDdRdth2.GetYaxis()->SetTitle("Delta t [ns]");
    IBDdRdth2.Draw("COLZ");

    //c3.SaveAs("/home/huangp/AntiNu/Plots/Reac_dRdt.png");

    
    
    // fitting each dr bin with gaussian
    // Define the fitting function (exponential decay)
    TF1* fit_func = new TF1("fit_func", "[0]*exp(-x/[1])", 0., 2000000.);
    TH2D* Analytical_IBDdRdth2 = new TH2D("Analytical_IBDdRdth2", "IBDdtdR", dRbin, 0., 2500, dtbin, 0., 2E6);
    TH2D* Acc_IBDdRdth2        = new TH2D("Acc_IBDdRdth2", "AccdtdR", dRbin, 0., 2500, dtbin, 0., 2E6);
    // Set parameter names and fix tau to 214150
    fit_func->SetParNames("Norm", "tau");
    fit_func->FixParameter(1, 214150);

    for( int ibinX =1; ibinX<= dRbin; ibinX++){
        double A  = IBDdRh1.GetBinContent(ibinX);
        double AA = Data_dRh1.GetBinContent(ibinX);
        FillIBD_2dhistocolumn(Analytical_IBDdRdth2,fit_func,A,ibinX);
        FillAcc_2dhistocolumn(Acc_IBDdRdth2,fit_func,AA,ibinX);
    }

    // Optionally, draw the result on a canvas
    TCanvas* c = new TCanvas("c", "Fit Slices Y", 800, 600);
    c->SetLogz();
    Analytical_IBDdRdth2->Scale(1. / Analytical_IBDdRdth2->Integral(), "nosw2");
    Analytical_IBDdRdth2->GetYaxis()->SetTitle("dt [ns]");
    Analytical_IBDdRdth2->GetXaxis()->SetTitle("dR [mm]");
    Analytical_IBDdRdth2->Draw("COLZ");
    c->SaveAs("finalIBDdtdr.png");

    TCanvas* cc = new TCanvas("cc", "Fit Slices ", 800, 600);
    cc->SetLogz();
    Acc_IBDdRdth2->GetYaxis()->SetTitle("dt [ns]");
    Acc_IBDdRdth2->GetXaxis()->SetTitle("dR [mm]");
    Acc_IBDdRdth2->Scale(1. / Acc_IBDdRdth2->Integral(), "nosw2");
    Acc_IBDdRdth2->Draw("COLZ");
    cc->SaveAs("finalACCdtdr.png");
    std::string output_root_address = "/home/huangp/AntiNu/Analysis/dtdrpdf.root";
    TFile* f_out = new TFile( output_root_address.c_str(),"RECREATE");
    Analytical_IBDdRdth2->Write();
    Acc_IBDdRdth2->Write();
    f_out->Close();
    delete f_out;

}