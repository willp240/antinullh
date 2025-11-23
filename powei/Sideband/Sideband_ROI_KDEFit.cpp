#include "RooKeysPdf.h"
#include "RooRealVar.h"
#include "RooWorkspace.h"
#include "RooDataSet.h"
#include "TFile.h"
#include "TCanvas.h"
#include "RooPlot.h"
#include "TRandom.h"
#include <iostream>

using namespace RooFit;

void Sideband_ROI_KDEFit() {
    // Load the workspace containing the RooKeysPdf
    TFile* file = TFile::Open("sidebandkdefit.root");
    RooWorkspace* w = (RooWorkspace*)file->Get("w_Po215");
    if (!w) {
        std::cerr << "Workspace not found in file!" << std::endl;
        return;
    }

    // Retrieve the RooKeysPdf and the observable
    RooKeysPdf* k1 = (RooKeysPdf*)w->pdf("k1");
    RooRealVar* delayedEcorr = (RooRealVar*)w->var("delayedEcorr");
    if (!k1 || !delayedEcorr) {
        std::cerr << "PDF or variable not found in workspace!" << std::endl;
        return;
    }

    // Define the range for sampling
    double rangeMin = 1.0, rangeMax = 2.5;
    delayedEcorr->setRange(rangeMin, rangeMax);

    // Generate a dataset by sampling from the KDE
    int nSamples = 100000; // Number of samples to generate
    RooDataSet* sampledData = k1->generate(*delayedEcorr, nSamples);
    RooKeysPdf newKDE("newKDE", "New KDE", *delayedEcorr, *sampledData, RooKeysPdf::NoMirror);
    
    RooDataSet scaled_data("scaled_data", "scaled dataset", RooArgSet(*delayedEcorr));

    Double_t Po214alphaPeak = 0.8432;Double_t Po215alphaPeak = 0.7910;
    Double_t scale_factor = Po214alphaPeak/Po215alphaPeak;
    for (int i = 0; i < sampledData->numEntries(); ++i) {
        const RooArgSet* entry = sampledData->get(i);
        double value = sampledData->get(i)->getRealValue("delayedEcorr");
        delayedEcorr->setVal(value * scale_factor); // Scale the value of x
        scaled_data.add(RooArgSet(*delayedEcorr));
    }
    RooKeysPdf newscaleKDE("newscaleKDE", "KDE for Po214", *delayedEcorr, scaled_data, RooKeysPdf::NoMirror);
    
    
    
    // Plot the sampled data and the original KDE
    RooPlot* frame = delayedEcorr->frame(Title("Sampling from KDE"));

    

    // Plot the sampled data
    sampledData->plotOn(frame, MarkerColor(kOrange), Name("Sampled Data"));
    scaled_data.plotOn(frame, MarkerColor(kRed), Name("Scaled Data"));
    // Plot the original PDF
    newKDE.plotOn(frame, LineColor(kGreen), Name("new KDE"));
    newscaleKDE.plotOn(frame, LineColor(kBlue), Name("New 214 Tail KDE"));
    //model->plotOn(frame, LineColor(kBlue), Name("Exponential PDF"));

    // Add legend
    TLegend* legend = new TLegend(0.6, 0.7, 0.9, 0.85);
    legend->SetTextSize(0.04);
    legend->SetBorderSize(0);
    legend->SetFillColorAlpha(0, 0); // Transparent background
    legend->AddEntry(frame->findObject("Sampled Data"), "Sampled Data", "P");
    legend->AddEntry(frame->findObject("Scaled Data"), "Scaled Data", "P");
    legend->AddEntry(frame->findObject("New KDE"), "215 Tail KDE", "L");
    legend->AddEntry(frame->findObject("New 214 Tail KDE"), "214 Tail KDE", "L");
    // Draw the plot
    TCanvas* c = new TCanvas("c", "Sampling from KDE", 800, 600);
    frame->Draw();
    legend->Draw();
    
    // Save the plot
    c->SaveAs("kde_sampling.png");


    // Optionally save the sampled dataset into a new workspace
    RooWorkspace *w2 = new RooWorkspace("w2", "workspace");
    w2->import(*sampledData);
    w2->import(scaled_data);
    //w2.import(*model);
    w2->import(newKDE);
    w2->import(newscaleKDE);
    w2->writeToFile("Sideband_ROI_KDEFit.root");
    w2->Print();
    // Cleanup
    delete sampledData;
    file->Close();
    delete file;
    delete c;

    std::cout << "Sampled data saved to sampledKDEWorkspace.root" << std::endl;
    std::cout << "Plot saved to kde_sampling.pdf" << std::endl;
    
}
