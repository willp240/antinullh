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
/// Script for plotting the unscaled PDFs inside a directory. The
/// inputs the directory name, and it loops through root files in that
/// directory and plots the PDF histograms to a PDF file, one PDF on
/// each page
///
/////////////////////////////////////////////////////////////////// */

void plot1DPDFs(std::string dirname)
{
    gStyle->SetOptStat(0);

    std::string outputfilename  = dirname + "/pdfplots.pdf";
    std::string outrootfilename = dirname + "/pdfplots.root";
    TFile *outfile = new TFile(outrootfilename.c_str(), "RECREATE");

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

    // Loop over files in the directory
    for (const auto &dirent : std::filesystem::directory_iterator(dirname))
    {
        const fs::path p = dirent.path();
        if (!dirent.is_regular_file()) continue;
        if (p.extension() != ".root")  continue;

        const std::string filePath = p.string();

        std::string histName = p.stem().string();

        // Base name before the last underscore (if any)
        std::string hname = histName;
        if (auto pos = histName.rfind('_'); pos != std::string::npos)
            hname.erase(pos);

        // Open ROOT file
        TFile *file = TFile::Open(filePath.c_str(), "READ");
        if (!file || file->IsZombie())
        {
            std::cerr << "Error opening file: " << filePath << std::endl;
            delete file;
            continue;
        }

        // Try to get histogram by full stem, then by base (before last '_')
        TH1D *h1 = nullptr;
        file->GetObject(histName.c_str(), h1);
        if (!h1 && hname != histName)
            file->GetObject(hname.c_str(), h1);

        if (!h1)
        {
            std::cerr << "No histogram named '" << histName
                      << "' or '" << hname
                      << "' found in file: " << filePath << std::endl;
            file->Close();
            delete file;
            continue;
        }

        // Draw and save page
        h1->SetLineWidth(2);
        h1->GetXaxis()->SetTitleSize(0.055);
        h1->GetYaxis()->SetTitleSize(0.055);
        h1->GetXaxis()->SetLabelSize(0.045);
        h1->GetYaxis()->SetLabelSize(0.045);
        h1->GetYaxis()->SetTitle("Probability");
        h1->SetTitle(h1->GetName());
        h1->Draw();
        gPad->SetGrid(1);
        gPad->Update();
        c1->Print(outputfilename.c_str());

        // Write the canvas into the output ROOT file
        outfile->cd();
        c1->Write(h1->GetName());

        // Clean up
        file->Close();
        delete file;
    }

    c1->Print((outputfilename + "]").c_str());
    c1->Close();
    outfile->Close();
    delete outfile;
}
