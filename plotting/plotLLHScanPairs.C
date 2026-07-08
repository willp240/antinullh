/* ///////////////////////////////////////////////////////////////////
///
/// Script for plotting LLH scans, with the same parameters for
/// different datasets plotted on the same canvas.
///
/// The user inputs the root file made by fixedOscLLHScan.
///
/// Histograms in the input file for the first dataset are looped
/// through and drawn. If a corresponding parameter exists for the
/// second dataset, the scan for it is drawn on the same canvas.
/// For parameters without a corresponding parameter in the other
/// dataset, the scan is drawn alone on a canvas.
///
/// The plots are drawn and saved in a pdf and root
/// file.
///
/////////////////////////////////////////////////////////////////// */

#include <TFile.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TKey.h>
#include <TClass.h>
#include <TString.h>
#include <iostream>
#include <map>
#include <cstring>
#include <cctype>
#include <string>

std::string trimTitle(const std::string &title)
{
    size_t first = 0;
    while (first < title.size() && std::isspace(static_cast<unsigned char>(title[first])))
        ++first;

    size_t last = title.size();
    while (last > first && std::isspace(static_cast<unsigned char>(title[last - 1])))
        --last;

    return title.substr(first, last - first);
}

std::string pairedAxisTitle(const char *title1, const char *title2)
{
    std::string a(title1 ? title1 : "");
    std::string b(title2 ? title2 : "");

    size_t prefix = 0;
    while (prefix < a.size() && prefix < b.size() && a[prefix] == b[prefix])
        ++prefix;

    size_t suffix = 0;
    while (suffix < a.size() - prefix && suffix < b.size() - prefix &&
           a[a.size() - 1 - suffix] == b[b.size() - 1 - suffix])
        ++suffix;

    std::string commonPrefix = trimTitle(a.substr(0, prefix));
    std::string commonSuffix = trimTitle(suffix > 0 ? a.substr(a.size() - suffix) : "");

    if (commonPrefix.empty())
        return a;
    if (commonSuffix.empty())
        return commonPrefix;
    return commonPrefix + " " + commonSuffix;
}

bool isOscillationScan(const TString &name)
{
    return name == "deltam21_full" || name == "sinsqtheta12_full";
}

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
    std::map<TString, int> cycles;
    TIter nextkey(f->GetListOfKeys());
    TKey *key;
    while ((key = (TKey *)nextkey()))
    {
        if (strcmp(key->GetClassName(), "TH1D") != 0)
            continue;
        TString histName = key->GetName();
        int cycle = key->GetCycle();
        if (cycles.find(histName) != cycles.end() && cycle < cycles.at(histName))
            continue;
        TH1D *h = (TH1D *)key->ReadObj();
        hists[histName] = h;
        cycles[histName] = cycle;
    }

    TString pdfName = outfiledir + "/llh_pair_plots.pdf";
    bool firstPage = true;

    // Loop through histograms
    for (const auto &kv : hists)
    {
        TString name = kv.first;
        std::cout << name << std::endl;

        // Skip bisMSB histos, handled via their PPO pair, and skip
        // correlated total/normalisation scans which have no dataset pair.
        if ((name.Contains("2_full") && !isOscillationScan(name)) || name.Contains("_norm_full"))
            continue;

        // Construct paired name
        TString pairName = name;
        pairName.ReplaceAll("_full", "2_full");

        TH1D *h1 = kv.second;
        TH1D *h2 = nullptr;

        bool hasPair = hists.find(pairName) != hists.end();
        if (!hasPair && !isOscillationScan(name))
            continue;

        if (hasPair)
            h2 = hists.at(pairName);

        // Make canvas
        TCanvas *c = new TCanvas(name, name, 1500, 1080);
        gStyle->SetOptStat(0);
        c->SetGrid();
        c->SetFrameLineWidth(2);
        c->cd();

        h1->SetLineColor(kBlue);
        h1->SetLineWidth(2);
        h1->SetTitle("");
        if (h2)
            h1->GetXaxis()->SetTitle(pairedAxisTitle(h1->GetXaxis()->GetTitle(), h2->GetXaxis()->GetTitle()).c_str());
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
        }
        if (h2)
            leg->Draw();

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
