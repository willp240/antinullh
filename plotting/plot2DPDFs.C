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
/// inputs are the directory name. It loops through root files in 
/// that directory and plots the PDF histograms to a PDF file, one
/// PDF on each page
///
/////////////////////////////////////////////////////////////////// */

void plot2DPDFs(std::string dirname)
{
    gStyle->SetOptStat(0);
    gStyle->SetPalette(51);

    std::vector<std::filesystem::path> pdfFiles;
    std::string outputfilename = dirname + "/pdfplots.pdf";
    std::string outrootfilename = dirname + "/pdfplots.root";
    TFile *outfile = new TFile(outrootfilename.c_str(), "RECREATE");

    // Add files from the pdfs directory
    for (const auto &entry : std::filesystem::directory_iterator(dirname.c_str()))
    {
        pdfFiles.push_back(entry.path());
    }

     // Print file name to first page
    TCanvas *c1 = new TCanvas("c1", "c1", 1500, 1080);
    c1->SetTopMargin(0.1);
    c1->SetBottomMargin(0.18);
    c1->SetLeftMargin(0.15);
    c1->SetRightMargin(0.13);
    c1->SetGrid();
    c1->SetFrameLineWidth(2);
    c1->Print((outputfilename + "[").c_str());
    c1->cd();

    // Loop over each file
    for (const auto &entry : pdfFiles)
    {

        std::string filePath = entry.string();
        std::filesystem::path pathObj(filePath);

        // Check if it's a .root file
        if (filePath.find(".root") == std::string::npos)
            continue;

        std::string histName = pathObj.filename().replace_extension("");

        TFile *file = TFile::Open(filePath.c_str(), "READ");
        if (!file || file->IsZombie())
        {
            std::cerr << "Error opening file: " << filePath << std::endl;
            delete file;
            continue;
        }

        // Get histogram
        TH2D *h1 = nullptr;
        file->GetObject(histName.c_str(), h1);

        if (!h1)
        {
            std::cerr << "No histogram " << histName << " found in file: " << filePath << std::endl;
            file->Close();
            delete file;
            continue;
        }

        h1->SetLineWidth(2);
        h1->GetXaxis()->SetTitleSize(0.055);
        h1->GetYaxis()->SetTitleSize(0.055);
        h1->GetXaxis()->SetLabelSize(0.045);
        h1->GetYaxis()->SetLabelSize(0.045);
        h1->SetTitle(histName.c_str());
        h1->Draw("colz");
        gPad->SetGrid(1);
        gPad->Update();
        c1->Print(outputfilename.c_str());

        outfile->cd();
        c1->Write(histName.c_str());
	gStyle->SetOptTitle(0);
	c1->SaveAs((dirname + histName + ".pdf").c_str());
	
    }

    c1->Print((outputfilename + "]").c_str());
    outfile->Close();

}
