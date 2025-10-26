// ROOT Headers
#include <TFile.h>
#include <TTree.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TStyle.h>

// c++ Headers
#include <sys/stat.h>

/* ///////////////////////////////////////////////////////////////////
///
/// Script for comparing the unscaled PDFs inside 2 directories. The
/// inputs the directory names, and it loops through root files in
/// the first directory and plots the PDF histogram. It checks the
/// same file and histogram exist in the second directory, and then
/// plots them on the same canvas. The plots are saved in a PDF file,
/// with one PDF comparison on each page
///
/////////////////////////////////////////////////////////////////// */

void compare1DPDFs(std::string dirname1, std::string dirname2, std::string label1, std::string label2)
{
    gStyle->SetOptStat(0);

    std::vector<std::filesystem::path> pdfFiles;
    std::string outputfilename = dirname1 + "/comppdfplots.pdf";
    std::string outrootfilename = dirname1 + "/comppdfplots.root";
    TFile *outfile = new TFile(outrootfilename.c_str(), "RECREATE");

    // Add files from the pdfs directory
    for (const auto &entry : std::filesystem::directory_iterator(dirname1.c_str()))
    {
        pdfFiles.push_back(entry.path());
    }

    // Print file name to first page
    TCanvas *c1 = new TCanvas("c1", "c1", 1500, 1080);
    c1->SetTopMargin(0.1);
    c1->SetBottomMargin(0.18);
    c1->SetLeftMargin(0.15);
    c1->SetRightMargin(0.09);
    c1->SetGrid();
    c1->SetFrameLineWidth(2);
    c1->Print((outputfilename + "[").c_str());
    c1->cd();

    // Loop over each file
    for (const auto &entry : pdfFiles)
    {

        std::string filePath1 = entry.string();
        std::filesystem::path pathObj(filePath1);
        // Check if it's a .root file
        if (filePath1.find(".root") == std::string::npos)
            continue;

        std::string histName = pathObj.filename().replace_extension("");
        TFile *file1 = TFile::Open(filePath1.c_str(), "READ");
        if (!file1 || file1->IsZombie())
        {
            std::cerr << "Error opening file: " << filePath1 << std::endl;
            delete file1;
            continue;
        }
        // Get histogram
        TH1D *h1 = nullptr;
        file1->GetObject(histName.c_str(), h1);
        if (!h1)
        {
            std::cerr << "No histogram " << histName << " found in file: " << filePath1 << std::endl;
            file1->Close();
            delete file1;
            continue;
        }

        std::string filePath2 = dirname2 + "/" + pathObj.filename().string();
        TFile *file2 = TFile::Open(filePath2.c_str(), "READ");
        if (!file2 || file2->IsZombie())
        {
            std::cerr << "Error opening file: " << filePath2 << std::endl;
            delete file2;
            continue;
        }
        // Get histogram
        TH1D *h2 = nullptr;
        file2->GetObject(histName.c_str(), h2);
        if (!h2)
        {
            std::cerr << "No histogram " << histName << " found in file: " << filePath2 << std::endl;
            file2->Close();
            delete file2;
            continue;
        }

        h1->SetMaximum(1.3 * h1->GetMaximum());
        h1->SetLineWidth(2);
        h1->SetLineColor(kBlue + 2);
        h1->GetXaxis()->SetTitleSize(0.055);
        h1->GetYaxis()->SetTitleSize(0.055);
        h1->GetXaxis()->SetLabelSize(0.045);
        h1->GetYaxis()->SetLabelSize(0.045);
        h1->GetYaxis()->SetTitle("Probability");
        h1->SetTitle("");
        h1->Draw();
        h2->Draw("same");
        h2->SetLineWidth(2);
        h2->SetLineStyle(2);
        h2->SetLineColor(kRed + 1);
        gPad->SetGrid(1);
        gPad->Update();

        TLegend *t1 = new TLegend(0.65, 0.73, 0.85, 0.88);
        t1->AddEntry(h1, label1.c_str(), "l");
        t1->AddEntry(h2, label2.c_str(), "l");
        t1->SetLineWidth(2);
        t1->SetTextFont(42);
        t1->Draw();

        c1->Print(outputfilename.c_str());
        outfile->cd();
        c1->Write(histName.c_str());
        c1->SaveAs((dirname1 + "/" + histName + ".pdf").c_str());
    }

    c1->Print((outputfilename + "]").c_str());
    outfile->Close();
}
