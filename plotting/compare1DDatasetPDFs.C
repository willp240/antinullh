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
/// Script for comparing the unscaled PDFs from 2 datasets. The
/// inputs are the directory name, dataset suffixes, and labels. It
/// loops through root files for the first dataset and checks the
/// same file and histogram exist for the second dataset, and then
/// plots them on the same canvas. The plots are saved in a PDF file,
/// with one PDF comparison on each page
///
/////////////////////////////////////////////////////////////////// */

void compare1DDatasetPDFs(std::string dirname, std::string suffix1, std::string suffix2, std::string label1, std::string label2)
{
    gStyle->SetOptStat(0);

    std::vector<std::filesystem::path> pdfFiles;
    std::string outputfilename = dirname + "/comppdfplots.pdf";
    std::string outrootfilename = dirname + "/comppdfplots.root";
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
        // Check if it's a .root file with first dataset suffix
        if (filePath1.find((suffix1+".root").c_str()) == std::string::npos)
            continue;

	std::string filename1 = pathObj.filename();
	std::string fullsuffix1 = suffix1 + ".root";
	std::string histname;

	size_t pos = filename1.rfind('_');  // find last underscore
	if (pos != std::string::npos) {
	  histname = filename1.substr(0, pos);   // keep only before the underscore
	}

        TFile *file1 = TFile::Open(filePath1.c_str(), "READ");
        if (!file1 || file1->IsZombie())
        {
            std::cerr << "Error opening file: " << filePath1 << std::endl;
            delete file1;
            continue;
        }
        // Get histogram
        TH1D *h1 = nullptr;
        file1->GetObject(histname.c_str(), h1);
        if (!h1)
        {
            std::cerr << "No histogram " << histname << " found in file: " << filePath1 << std::endl;
            file1->Close();
            delete file1;
            continue;
        }

	std::string histname2 = histname + "2";
        std::string filePath2 = dirname + "/" + histname2 + "_" + suffix2 + ".root";
        TFile *file2 = TFile::Open(filePath2.c_str(), "READ");
        if (!file2 || file2->IsZombie())
        {
            std::cerr << "Error opening file: " << filePath2 << std::endl;
            delete file2;
            continue;
        }
        // Get histogram
        TH1D *h2 = nullptr;
        file2->GetObject(histname2.c_str(), h2);
        if (!h2)
        {
            std::cerr << "No histogram " << histname2 << " found in file: " << filePath2 << std::endl;
            file2->Close();
            delete file2;
            continue;
        }

        h1->GetYaxis()->SetRangeUser(0, (1.3 * h1->GetMaximum()));
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
        c1->Write(histname.c_str());
        c1->SaveAs((dirname + "/" + histname + ".pdf").c_str());
    }

    c1->Print((outputfilename + "]").c_str());
    outfile->Close();
}
