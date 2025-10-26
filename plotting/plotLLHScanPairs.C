#include <TFile.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TKey.h>
#include <TClass.h>
#include <TString.h>
#include <iostream>
#include <map>

void plotLLHScanPairs(const char *filename = "llh_scan.root")
{

    // Open file
    std::filesystem::path filepath(filename);
    std::string outputfilename = filepath.replace_extension("plots.pdf").string();
    std::string outrootfilename = filepath.replace_extension(".root").string();
    std::string outfiledir = filepath.parent_path().string();
    TFile *outfile = new TFile(outrootfilename.c_str(), "RECREATE");

    // Open file
    TFile *f = TFile::Open(filename);
    if (!f || f->IsZombie())
    {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }

    // Collect all TH1D histograms
    std::map<TString, TH1D *> hists;
    TIter nextkey(f->GetListOfKeys());
    TKey *key;
    while ((key = (TKey *)nextkey()))
    {
        if (strcmp(key->GetClassName(), "TH1D") != 0)
            continue;
        TH1D *h = (TH1D *)key->ReadObj();
        hists[h->GetName()] = h;
    }

    TString pdfName = outfiledir + "/llh_pair_plots.pdf";
    bool firstPage = true;

    // Loop through histograms
    for (const auto &kv : hists)
    {
        TString name = kv.first;
        std::cout << name << std::endl;

        // Skip "_2_full" histos, handle via their pair
        if (name.Contains("2_full") && !name.Contains("theta"))
            continue;

        // Construct paired name
        TString pairName = name;
        pairName.ReplaceAll("_full", "2_full");

        TH1D *h1 = kv.second;
        TH1D *h2 = nullptr;

        if (hists.find(pairName) != hists.end())
        {
            h2 = hists.at(pairName);
        }

        // Make canvas
        TCanvas *c = new TCanvas(name, name, 1500, 1080);
        gStyle->SetOptStat(0);
        c->SetGrid();
        c->SetFrameLineWidth(2);
        c->cd();

        h1->SetLineColor(kBlue);
        h1->SetLineWidth(2);
        h1->SetTitle("");
        h1->Draw("HIST");

        TLegend *leg = new TLegend(0.4, 0.6, 0.6, 0.8);
        leg->SetLineWidth(2);
        leg->AddEntry(h1, "PPO", "l");

        if (h2)
        {
            h2->SetLineColor(kRed);
            h2->SetLineStyle(2);
            h2->SetLineWidth(2);
            h2->Draw("HIST SAME");
            leg->AddEntry(h2, "bisMSB", "l");
            leg->Draw();
        }

        // Save individual pdf
        TString pngName = outfiledir + "/" + name + ".pdf";
        c->SaveAs(pngName);

        // Append to PDF (multi-page)
        if (firstPage)
        {
            c->SaveAs(pdfName + "("); // open
            firstPage = false;
        }
        else
        {
            c->SaveAs(pdfName);
        }

        outfile->cd();
        c->Write();

        delete c;
    }

    // Close PDF
    TCanvas *dummy = new TCanvas("dummy", "", 800, 600);
    dummy->SaveAs(pdfName + ")");

    std::cout << "Saved multi-page PDF: " << pdfName << std::endl;
    f->Close();
}
