// Makes prefit-postfit comparison plots and correlations
// Loops over each parameter and plots a 1D distribution of
// parameter value at each step. Calculates mean and error in
// 4 different ways, and also plots the asimov value. Then
// makes a plot of means and errors for all parameters. Also
// the pre and postfit MC distributions, projections of these,
// and ratios between them
// ./make_plots will show you the options
#include <iostream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <TMath.h>
#include <string>

#include "TObjArray.h"
#include "TChain.h"
#include "TFile.h"
#include "TBranch.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TLegend.h"
#include "TString.h"
#include "TStyle.h"
#include "TH1.h"
#include "TF1.h"
#include "TH2.h"
#include "TGraphErrors.h"
#include "TVectorD.h"
#include "TColor.h"
#include "TKey.h"
#include "TPaveText.h"

#include <FitConfigLoader.hh>
#include <FitConfig.hh>

// Should put these in header file
std::vector<std::string> *parNames;
std::vector<double> *parRates;
std::vector<TString> bnames;
typedef std::map<std::string, double> ParMap;
ParMap constrMeans;
ParMap constrSigmas;
ParMap mins;
ParMap maxs;
ParMap asimovRates;
std::string mcmcConfigFile;
std::string asmvDistsFileName;
std::string scaledPostfitDistFileName;

TH1D *MakePrefit(int npar);

// Get the highest posterior density value and errors
void GetHPD(TH1D *const post, double &central, double &error, double &error_pos, double &error_neg);
// Get arithmetic mean and sigma of histogram
void GetArithmetic(TH1D *const hpost, double &mean, double &error);
// Fit gaussian and get mean and sigma
void GetGaussian(TH1D *&hpost, TF1 *&gauss, double &central, double &error);
// Get parameter value at max LLH step
double GetMaxLLHPar(TTree *chain, int maxStep, std::string parname);
// Get max LLH step
int GetMaxLLHStep(TTree *chain);

// Main function with loop over all params
void make_plots(std::string inputfile, bool correlations = false, bool traces = true);

// Get bounds for each parameter
void GetParLimits(std::string paramName, double &central, double &prior, double &down_error, double &up_error);

// Load up values from input files
void LoadInputVals();
// a nice palette for ratio plots
void BlueRedPalette();

// Make 1D prefit/postfit comparisons
TCanvas *CompareProjections(TH1D *prefit, TH1D *postfit);

// How much burn in should we take off? As a fraction of total steps i.e BurnInCut=5 means reject the first 1/5th of steps
int BurnInCut = 5;

using namespace bbfit;

int main(int argc, char *argv[])
{

  if (argc != 5 && argc != 6 && argc != 7)
  {
    std::cerr << "./make_plots root_file_to_analyse.root asimovDistFile scaledPostfitDistFile mcmcConfigFile correlations traces" << std::endl;
    exit(-1);
  }

  std::string filename = argv[1];
  asmvDistsFileName = argv[2];
  scaledPostfitDistFileName = argv[3];
  mcmcConfigFile = argv[4];

  // Plotting correlations gives lots of plots and takes a bit of time, so only do if user wants to
  bool correlations = false;
  if (argc == 6)
  {
    correlations = true;
  }
  // Plotting traces gives lots of plots and takes a bit of time, so only do if user wants to
  bool traces = false;
  if (argc == 7)
  {
    traces = true;
  }

  make_plots(filename, correlations, traces);

  return 0;
}

void make_plots(std::string inputFile, bool correlations, bool traces)
{

  // do we want to draw correlations and traces or not
  std::cout << "File for study:       " << inputFile << std::endl;
  std::cout << "Draw correlations?    " << correlations << std::endl;
  std::cout << "Draw traces?           " << traces << std::endl;

  // Load up prefit values
  LoadInputVals();

  // Open the chain
  TChain *chain = new TChain("posteriors", "");
  chain->Add(inputFile.c_str());
  // MCMC Burn-in cut
  int cut = chain->GetMaximum("Step") / BurnInCut;
  std::stringstream ss;
  ss << "Step > " << cut;

  std::cout << "MCMC has " << chain->GetMaximum("Step") << " steps, burn-in is set to " << cut << std::endl;
  std::string stepcut = ss.str();

  // Get the list of branches
  TObjArray *brlis = (TObjArray *)chain->GetListOfBranches();
  // Get the number of branches
  int nbr = brlis->GetEntries();
  std::cout << "# of branches: " << nbr << std::endl;

  // Have a counter for how many parameters we have
  int npar = 0;

  int maxLLHStep = GetMaxLLHStep(chain);
  std::vector<bool> systFlag;

  // Loop over the number of branches
  chain->SetBranchStatus("*", false);
  chain->SetBranchStatus("Step", true);
  for (int i = 0; i < nbr; i++)
  {

    // Get the TBranch and its name
    TBranch *br = (TBranch *)brlis->At(i);
    TString bname = br->GetName();

    if (bname.BeginsWith("LogL"))
      continue;
    if (bname.BeginsWith("Accepted"))
      continue;
    if (bname.BeginsWith("Step"))
      continue;
    if (bname.BeginsWith("energy_") || bname.BeginsWith("r_"))
      systFlag.push_back(true);
    else
      systFlag.push_back(false);

    chain->SetBranchStatus(bname, true);
    bnames.push_back(bname);
    npar++;
  }
  std::cout << "# of parameters: " << npar << std::endl;

  // Get first entry in chain
  chain->GetEntry(0);

  gStyle->SetOptFit(111);

  // Open a TCanvas to write the posterior onto
  TCanvas *c0 = new TCanvas("c0", "c0", 0, 0, 1600, 1024);
  c0->SetGrid();
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetNumberContours(255);
  c0->SetTickx();
  c0->SetTicky();
  c0->SetBottomMargin(0.1);
  c0->SetTopMargin(0.02);
  c0->SetRightMargin(0.1);
  c0->SetLeftMargin(0.10);

  // Make sure we can read files located anywhere and strip the .root ending
  inputFile = inputFile.substr(0, inputFile.find(".root"));

  TString canvasname = inputFile;
  // Append if we're drawing correlations
  if (correlations)
  {
    canvasname += "_plotCorrelations.pdf[";
  }
  else
  {
    canvasname += "_plotParameters.pdf[";
  }
  c0->Print(canvasname);

  // Once the pdf file is open no longer need to bracket
  canvasname.ReplaceAll("[", "");

  TF1 *gauss = new TF1("gauss", "[0]/sqrt(2.0*TMath::Pi())/[2]*TMath::Exp(-0.5*pow(x-[1],2)/[2]/[2])", -5, 5);
  gauss->SetLineWidth(2);
  gauss->SetLineColor(kGreen + 3);

  // Some TVectors with the means, errors and Gaussian parameters of the PDFs
  TVectorD *mean_vec = new TVectorD(npar);
  TVectorD *err_vec = new TVectorD(npar);
  TVectorD *gaus_mean_vec = new TVectorD(npar);
  TVectorD *gaus_err_vec = new TVectorD(npar);
  TVectorD *HPD_mean_vec = new TVectorD(npar);
  TVectorD *HPD_err_p_vec = new TVectorD(npar);
  TVectorD *HPD_err_vec = new TVectorD(npar);
  TVectorD *HPD_err_m_vec = new TVectorD(npar);
  TVectorD *HLLH_mean_vec = new TVectorD(npar);

  TMatrixT<double> *correlation = new TMatrixT<double>(npar, npar);
  for (int i = 0; i < npar; ++i)
  {
    for (int j = 0; j < npar; ++j)
    {
      (*correlation)(i, j) = 0.;
    }
  }

  // Output file to write to
  TString rootfilename = inputFile;
  if (correlations)
  {
    rootfilename += "_plotCorrelations.root";
  }
  else
  {
    rootfilename += "_plotParameters.root";
  }

  // The output file
  TFile *file = new TFile(rootfilename, "RECREATE");
  file->cd();
  // Make a directory with the parameters
  TDirectory *params = file->mkdir("params");
  TDirectory *corr = gDirectory;
  if (correlations)
  {
    corr = file->mkdir("corr");
  }

  // Now loop over all parameters
  for (int i = 0; i < npar; ++i)
  {

    file->cd();
    // Number of bins in 1D plot
    int nbins = 70;
    double asimovLine = 0.0;

    // Get bounds and asimov rates for parameter
    std::string tempString = std::string(bnames.at(i));
    double central, prior, down, up;
    GetParLimits(tempString, central, prior, down, up);
    asimovLine = central;

    file->cd();

    // Get the maximum and minimum for the parameter
    chain->Draw(bnames.at(i), stepcut.c_str());
    TH1 *htemp = (TH1 *)gPad->GetPrimitive("htemp");
    double maximum = htemp->GetXaxis()->GetBinUpEdge(htemp->GetXaxis()->GetNbins());
    double minimum = htemp->GetXaxis()->GetBinLowEdge(1);

    if (minimum < 0 && !systFlag.at(i))
    {
      htemp->GetXaxis()->SetRangeUser(0, maximum);
    }
    // This holds the posterior density
    TH1D *hpost = new TH1D(bnames.at(i), bnames.at(i), nbins, minimum, maximum);
    hpost->SetMinimum(0);
    hpost->GetYaxis()->SetTitle("Steps");
    hpost->GetYaxis()->SetNoExponent(false);
    hpost->SetTitle(bnames.at(i));

    // Project bnames.at(i) onto hpost, applying stepcut
    chain->Project(bnames.at(i), bnames.at(i), stepcut.c_str());

    // Apply one smoothing
    hpost->Smooth();

    // Get the characteristics of the hpost
    double mean, rms;
    GetArithmetic(hpost, mean, rms);
    double peakval, sigma_p, sigma_m, sigma_hpd;
    GetHPD(hpost, peakval, sigma_hpd, sigma_p, sigma_m);
    double gauss_mean, gauss_rms;
    GetGaussian(hpost, gauss, gauss_mean, gauss_rms);
    double maxllh = GetMaxLLHPar(chain, maxLLHStep, bnames.at(i).Data());
    double maxllh_norm;

    std::cout << i << ": " << mean << " +/- " << rms << " (" << peakval << "+/-" << sigma_hpd << " + " << sigma_p << " - " << sigma_m << ")" << " (" << gauss_mean << "+/-" << gauss_rms << ")" << std::endl;

    TLine *hpd = new TLine(peakval, hpost->GetMinimum(), peakval, hpost->GetMaximum());
    hpd->SetLineColor(kRed);
    hpd->SetLineWidth(2);
    hpd->SetLineStyle(kSolid);

    // Make the legend
    TLegend *leg = new TLegend(0.12, 0.6, 0.6, 0.97);
    leg->SetTextSize(0.04);
    leg->AddEntry(hpost, Form("#splitline{PDF}{#mu = %.2f, #sigma = %.2f}", hpost->GetMean(), hpost->GetRMS()), "l");
    leg->AddEntry(gauss, Form("#splitline{Gauss}{#mu = %.2f, #sigma = %.2f}", gauss->GetParameter(1), gauss->GetParameter(2)), "l");
    leg->AddEntry(hpd, Form("#splitline{HPD}{#mu = %.2f, #sigma = %.2f (+%.2f-%.2f)}", peakval, sigma_hpd, sigma_p, sigma_m), "l");

    // Now take values as a ratio to asimov rate. If asimov rate is 0, set central value to 1 for plotting
    if (central != 0)
    {
      mean = mean / central;
      rms = rms / central;
      gauss_mean = gauss_mean / central;
      gauss_rms = gauss_rms / central;
      peakval = peakval / central;
      sigma_hpd = sigma_hpd / central;
      sigma_p = sigma_p / central;
      sigma_m = sigma_m / central;
      maxllh_norm = maxllh / central;
    }
    else if (central == 0)
    {
      mean = mean + 1.0;
      gauss_mean = gauss_mean + 1.0;
      peakval = peakval + 1.0;
      maxllh_norm = maxllh;
    }

    (*mean_vec)(i) = mean;
    (*err_vec)(i) = rms;
    (*gaus_mean_vec)(i) = gauss_mean;
    (*gaus_err_vec)(i) = gauss_rms;
    (*HPD_mean_vec)(i) = peakval;
    (*HPD_err_p_vec)(i) = sigma_p;
    (*HPD_err_vec)(i) = sigma_hpd;
    (*HPD_err_m_vec)(i) = sigma_m;
    (*correlation)(i, i) = 1.0;
    (*HLLH_mean_vec)(i) = maxllh_norm;

    hpost->SetLineWidth(2);
    hpost->SetLineColor(kBlack);
    hpost->SetMaximum(hpost->GetMaximum() * 1.5);
    hpost->SetTitle(tempString.c_str());
    if (systFlag.at(i))
      hpost->GetXaxis()->SetTitle(("Value of " + std::string(hpost->GetTitle())).c_str());
    else
      hpost->GetXaxis()->SetTitle(("Number of " + std::string(hpost->GetTitle()) + " Events").c_str());

    // Now make the TLine for the asimov
    TLine *asimov = new TLine(asimovLine, hpost->GetMinimum(), asimovLine, hpost->GetMaximum());
    asimov->SetLineColor(kBlue);
    asimov->SetLineWidth(2);
    asimov->SetLineStyle(kDashed);

    // And the line for parameter value in the single step that has biggest LLH
    TLine *maxllhline = new TLine(maxllh, hpost->GetMinimum(), maxllh, hpost->GetMaximum());
    maxllhline->SetLineColor(kMagenta);
    maxllhline->SetLineWidth(2);
    maxllhline->SetLineStyle(kDashed);

    // And draw it all
    hpost->Draw();
    hpd->Draw("same");
    asimov->Draw("same");
    maxllhline->Draw("same");

    if (prior)
      leg->AddEntry(asimov, Form("#splitline{Asimov}{x = %.2f}", asimovLine), "l");
    else
      leg->AddEntry(asimov, Form("#splitline{Asimov}{x = %.2f}", asimovLine), "l");
    leg->AddEntry(maxllhline, Form("#splitline{Max LLH}{x = %.2f}", maxllh), "l");
    leg->SetLineColor(0);
    leg->SetLineStyle(0);
    leg->SetFillColor(0);
    leg->SetFillStyle(0);
    leg->Draw("same");

    // Write to file
    c0->SetName(hpost->GetName());
    c0->SetTitle(hpost->GetTitle());
    c0->Print(canvasname);

    // cd into params directory in root file
    file->cd();
    params->cd();
    c0->Write();

    if (traces)
    {
      TH2D *htrace = new TH2D(bnames.at(i) + "_trace", bnames.at(i) + "_trace", 100, 0, chain->GetMaximum("Step"), 100, minimum, maximum);
      htrace->GetXaxis()->SetTitle("Step");
      htrace->GetYaxis()->SetTitle(bnames.at(i));
      htrace->SetTitle(bnames.at(i) + "_trace");

      // Project bnames.at(i) onto htrace, applying stepcut
      chain->Project(bnames.at(i) + "_trace", bnames.at(i) + ":Step");

      gStyle->SetPalette(51);
      htrace->Draw("colz");
      c0->SetName(htrace->GetName());
      c0->SetTitle(htrace->GetTitle());
      c0->Print(canvasname);
    }

    // If we're interested in drawing the correlations need to invoke another for loop
    // Can surely improve this with less repeated code but it will do for now
    if (correlations)
    {

      // Loop over the other parameters to get the correlations
      for (int j = 0; j <= i; j++)
      {

        // Skip the diagonal elements which we've already done above
        if (j == i)
          continue;

        std::string tempString = std::string(bnames[j]);
        // Get the parameter name
        double central, prior, down, up;
        GetParLimits(tempString, central, prior, down, up);
        asimovLine = central;

        file->cd();

        int nbins = 70;

        TString drawcmd = bnames.at(j) + ":" + bnames.at(i);

        chain->Draw(bnames.at(j), stepcut.c_str());
        TH1 *htemp2 = (TH1 *)gPad->GetPrimitive("htemp");
        double maximum2 = htemp2->GetXaxis()->GetBinUpEdge(htemp2->GetXaxis()->GetNbins());
        double minimum2 = htemp2->GetXaxis()->GetBinLowEdge(1);

        // TH2F to hold the correlation
        TH2F *hpost2 = new TH2F(drawcmd, drawcmd, hpost->GetNbinsX(), hpost->GetBinLowEdge(0), hpost->GetBinLowEdge(hpost->GetNbinsX() + 1), nbins, minimum2, maximum2);
        hpost2->SetMinimum(0);
        hpost2->GetXaxis()->SetTitle(hpost->GetXaxis()->GetTitle());
        hpost2->GetYaxis()->SetTitle(("Number of " + std::string(bnames.at(j)) + " Events").c_str());
        hpost2->GetYaxis()->SetTitleOffset(1.2);
        hpost2->GetZaxis()->SetTitle("Steps");
        std::string plotTitle = hpost->GetXaxis()->GetTitle();
        plotTitle += " vs " + tempString;
        hpost2->SetTitle(plotTitle.c_str());

        // The draw command we want, i.e. draw param j vs param i
        chain->Project(drawcmd, drawcmd, stepcut.c_str());
        hpost2->Draw("colz");

        (*correlation)(i, j) = hpost2->GetCorrelationFactor();
        (*correlation)(j, i) = (*correlation)(i, j);
        std::cout << std::setw(10) << "correlation:" << (*correlation)(i, j) << std::endl;

        gStyle->SetPalette(51);
        c0->SetName(hpost2->GetName());
        c0->SetTitle(hpost2->GetTitle());
        c0->Print(canvasname);

        // Write it to root file
        file->cd();
        corr->cd();
        hpost2->Write();

        delete hpost2;

      } // End for (j = 0; j <= i; ++j)
    }
    delete hpost;
    // delete htrace;
    delete asimov;
    delete hpd;
    delete leg;

  } // End for npar

  // Make the prefit plot
  TH1D *prefit = MakePrefit(npar);

  // cd into the output file
  file->cd();

  // Make a TH1D of the central values and the errors
  TH1D *paramPlot = new TH1D("paramPlot", "paramPlot", npar, 0, npar);
  paramPlot->SetName("postfitparams");
  paramPlot->SetTitle(stepcut.c_str());
  paramPlot->SetFillColor(kBlack);
  paramPlot->SetMarkerColor(paramPlot->GetFillColor());
  paramPlot->SetMarkerStyle(108);
  paramPlot->SetLineColor(paramPlot->GetFillColor());
  paramPlot->SetLineWidth(3);
  paramPlot->SetMarkerSize(prefit->GetMarkerSize());

  // Same but with Gaussian output
  TH1D *paramPlot_gauss = (TH1D *)(paramPlot->Clone());
  paramPlot_gauss->SetMarkerColor(kGreen + 3);
  paramPlot_gauss->SetMarkerStyle(23);
  paramPlot_gauss->SetLineWidth(2);
  paramPlot_gauss->SetMarkerSize((prefit->GetMarkerSize()) * 0.75);
  paramPlot_gauss->SetFillColor(paramPlot_gauss->GetMarkerColor());
  paramPlot_gauss->SetFillStyle(3244);
  paramPlot_gauss->SetLineColor(paramPlot_gauss->GetMarkerColor());

  // Same but with HPD output
  TH1D *paramPlot_HPD = (TH1D *)(paramPlot->Clone());
  paramPlot_HPD->SetMarkerColor(kRed);
  paramPlot_HPD->SetMarkerStyle(25);
  paramPlot_HPD->SetLineWidth(2);
  paramPlot_HPD->SetMarkerSize((prefit->GetMarkerSize()) * 0.5);
  paramPlot_HPD->SetFillColor(kRed);
  paramPlot_HPD->SetFillStyle(3154);
  paramPlot_HPD->SetLineColor(paramPlot_HPD->GetMarkerColor());

  // And for highest LLH
  TH1D *paramPlot_HLLH = (TH1D *)(paramPlot->Clone());
  paramPlot_HLLH->SetFillColor(kMagenta);
  paramPlot_HLLH->SetMarkerColor(paramPlot_HLLH->GetFillColor());
  paramPlot_HLLH->SetMarkerStyle(108);
  paramPlot_HLLH->SetLineColor(paramPlot_HLLH->GetFillColor());
  paramPlot_HLLH->SetLineWidth(3);
  paramPlot_HLLH->SetMarkerSize(prefit->GetMarkerSize());

  // Set labels and data
  for (int i = 0; i < npar; ++i)
  {
    paramPlot->SetBinContent(i + 1, (*mean_vec)(i));
    paramPlot->SetBinError(i + 1, (*err_vec)(i));

    paramPlot_gauss->SetBinContent(i + 1, (*gaus_mean_vec)(i));
    paramPlot_gauss->SetBinError(i + 1, (*gaus_err_vec)(i));

    paramPlot_HPD->SetBinContent(i + 1, (*HPD_mean_vec)(i));
    double error = (*HPD_err_vec)(i);
    paramPlot_HPD->SetBinError(i + 1, error);

    paramPlot_HLLH->SetBinContent(i + 1, (*HLLH_mean_vec)(i));
    // This is a bit of a fudge. Don't have error for HLLH as it's just a parameter value at a single step. But if the error is set to 0 the plot doesn't show up (not sure why, prolly a root thing). Could set to a small but finite number, say 0.01, but as we plot over log y scale the line draws thicker for lower parameter values. So to show a 0 error (to the eye) and have constant line size, we set the error to be a very small constant fraction of the bin content
    paramPlot_HLLH->SetBinError(i + 1, 0.01 * paramPlot_HLLH->GetBinContent(i + 1));
  }

  // Make a TLegend
  TLegend *CompLeg = new TLegend(0.12, 0.80, 0.5, 0.95);
  CompLeg->AddEntry(prefit, "Prior", "fp");
  CompLeg->AddEntry(paramPlot_HPD, "Postfit HPD", "fp");
  CompLeg->AddEntry(paramPlot, "Postfit PDF", "lep");
  CompLeg->AddEntry(paramPlot_gauss, "Postfit Gauss", "fp");
  CompLeg->AddEntry(paramPlot_HLLH, "Postfit HLLH", "lep");
  CompLeg->SetFillColor(0);
  CompLeg->SetFillStyle(0);
  CompLeg->SetLineWidth(0);
  CompLeg->SetLineStyle(0);
  CompLeg->SetBorderSize(0);

  file->cd();
  c0->SetBottomMargin(0.2);

  file->cd();
  prefit->GetYaxis()->SetRangeUser(0.1, 10000.0);
  gPad->SetLogy();
  // TODO: config on log axis
  // prefit->GetYaxis()->SetRangeUser(0.001, 2.0);
  prefit->GetXaxis()->SetTitle("");
  prefit->GetXaxis()->SetRangeUser(0, npar);
  prefit->GetXaxis()->LabelsOption("v");

  // Plot and write all the full parameter histograms
  paramPlot->GetXaxis()->SetRangeUser(0, npar);
  paramPlot_gauss->GetXaxis()->SetRangeUser(0, npar);
  paramPlot_HPD->GetXaxis()->SetRangeUser(0, npar);
  paramPlot_HLLH->GetXaxis()->SetRangeUser(0, npar);

  prefit->Write("param_prefit");
  paramPlot->Write("param");
  paramPlot_gauss->Write("param_gaus");
  paramPlot_HPD->Write("param_HPD");
  paramPlot_HPD->Write("param_HLLH");

  prefit->Draw("e2");
  paramPlot_gauss->Draw("e2, same");
  paramPlot_HPD->Draw("e2, same");
  paramPlot->Draw("same");
  paramPlot_HLLH->Draw("e2, same");

  CompLeg->SetX1NDC(0.33);
  CompLeg->SetX2NDC(0.80);
  CompLeg->SetY1NDC(0.20);
  CompLeg->SetY2NDC(0.50);

  CompLeg->Draw("same");
  c0->Write("param_canv");
  c0->Print(canvasname);
  c0->Clear();

  delete CompLeg;

  c0->SetLeftMargin(0.1);
  c0->SetBottomMargin(0.1);

  // Now do printing for correlations
  TH2D *hCorr = NULL;
  if (correlations)
  {
    gPad->SetLogy(0);
    // The correlation
    hCorr = new TH2D("hCorr", "hCorr", npar, 0, npar, npar, 0, npar);
    hCorr->GetZaxis()->SetTitle("Correlation");
    hCorr->SetMinimum(-1);
    hCorr->SetMaximum(1);
    hCorr->GetXaxis()->SetLabelSize(0.020);
    hCorr->GetYaxis()->SetLabelSize(0.020);

    // Loop over the correlation matrix entries
    for (int i = 0; i < npar; i++)
    {

      hCorr->GetXaxis()->SetBinLabel(i + 1, prefit->GetXaxis()->GetBinLabel(i + 1));

      for (int j = 0; j < npar; j++)
      {

        hCorr->GetYaxis()->SetBinLabel(j + 1, prefit->GetXaxis()->GetBinLabel(j + 1));

        // The value of the correlation
        double corr = (*correlation)(i, j);

        hCorr->SetBinContent(i + 1, j + 1, corr);
      }
    }

    BlueRedPalette();

    c0->cd();
    c0->Clear();
    hCorr->Draw("colz");
    c0->SetRightMargin(0.15);
    c0->Print(canvasname);

    file->cd();
    hCorr->Write("postfit_corr_plot");
  }

  // Read prefit event distribution
  TFile *asmvDistsFile = new TFile(asmvDistsFileName.c_str(), "READ");
  asmvDistsFile->cd();
  TH2D *prefitDist;
  TIter next(asmvDistsFile->GetListOfKeys());
  TKey *key;
  while ((key = (TKey *)next()))
  {
    std::string classname = std::string(key->GetClassName());
    if (classname == "TH2D")
      prefitDist = (TH2D *)key->ReadObj();
  }

  // Read postfit event distribution
  TFile *scaledPostfitDistFile = new TFile(scaledPostfitDistFileName.c_str(), "READ");
  scaledPostfitDistFile->cd();
  TH2D *postfitDist;
  TIter next2(scaledPostfitDistFile->GetListOfKeys());
  TKey *key2;

  while ((key2 = (TKey *)next2()))
  {
    std::string classname = std::string(key2->GetClassName());
    if (classname == "TH2D")
      postfitDist = (TH2D *)key2->ReadObj();
  }

  gStyle->SetPalette(51);
  gStyle->SetOptTitle(1);
  gPad->SetLogy(0);
  c0->SetTopMargin(0.1);
  c0->SetRightMargin(0.13);

  // Draw prefit and postfit event distributions
  prefitDist->SetTitle("Prefit Distribution");
  prefitDist->Draw("colz");
  c0->Print(canvasname);

  postfitDist->SetTitle("Postfit Distribution");
  postfitDist->Draw("colz");
  c0->Print(canvasname);

  BlueRedPalette();

  // Take ratio prefit over postfit
  TH2D *ratioplot = (TH2D *)prefitDist->Clone();
  ratioplot->SetTitle("Prefit/Postfit");
  ratioplot->Divide(postfitDist);
  ratioplot->GetZaxis()->SetRangeUser(0.7, 1.3);
  ratioplot->Draw("colz");
  c0->Print(canvasname);

  // And let's look at the the absolute difference as well
  TH2D *diffplot = (TH2D *)prefitDist->Clone();
  diffplot->SetTitle("Prefit - Postfit");
  diffplot->Add(postfitDist, -1);
  // A bit fiddly here, but basically want z axis to be symmetric around 0
  double maxZ = diffplot->GetMaximum();
  double minZ = diffplot->GetMinimum();
  double absMaxZ;
  if (maxZ > -minZ)
    absMaxZ = maxZ;
  else
    absMaxZ = -minZ;
  diffplot->GetZaxis()->SetRangeUser(-absMaxZ, absMaxZ);
  diffplot->Draw("colz");
  c0->Print(canvasname);

  // Now make 1D projections
  TH1D *PrefitE = (TH1D *)prefitDist->ProjectionX()->Clone();
  PrefitE->SetName("Prefit_e");
  TH1D *PrefitR = (TH1D *)prefitDist->ProjectionY()->Clone();
  PrefitE->SetName("Prefit_r");
  TH1D *PostfitE = (TH1D *)postfitDist->ProjectionX()->Clone();
  PrefitE->SetName("Postfit_e");
  TH1D *PostfitR = (TH1D *)postfitDist->ProjectionY()->Clone();
  PrefitE->SetName("Postfit_r");

  // Make comparisons of these 1D projections
  TCanvas *e1Ds = CompareProjections(PrefitE, PostfitE);
  e1Ds->SetName("e1Ds");
  e1Ds->SetTitle("e1Ds");
  TCanvas *r1Ds = CompareProjections(PrefitR, PostfitR);
  r1Ds->SetName("r1Ds");
  r1Ds->SetTitle("r1Ds");

  e1Ds->Print(canvasname);
  r1Ds->Print(canvasname);

  // Then close the pdf file
  std::cout << "Closing pdf " << canvasname << std::endl;
  canvasname += "]";
  c0->Print(canvasname);

  // Write all the nice vectors
  file->cd();
  mean_vec->Write("postfit_params_arit");
  err_vec->Write("postfit_errors_arit");
  gaus_mean_vec->Write("postfit_params_gauss");
  gaus_err_vec->Write("postfit_errors_gauss");
  HPD_mean_vec->Write("postfit_params_HPD");
  HPD_err_vec->Write("postfit_errors_HPD");
  HPD_err_p_vec->Write("postfit_errors_HPD_pos");
  HPD_err_m_vec->Write("postfit_errors_HPD_neg");
  HLLH_mean_vec->Write("postfit_params_HLLH");
  correlation->Write("postfit_corr");
  postfitDist->Write("postfit_dist");
  prefitDist->Write("prefit_dist");
  ratioplot->Write("prefit_postfit_ratio");
  diffplot->Write("prefit_postfit_diff");
  e1Ds->Write("e_comparisons");
  r1Ds->Write("r_comparisons");

  file->Close();
}

// **************************
// Get the highest posterior density from a TH1D
void GetHPD(TH1D *const hpost, double &central, double &error, double &error_pos, double &error_neg)
{
  // **************************

  // Get the bin which has the largest posterior density
  int MaxBin = hpost->GetMaximumBin();
  // And it's value
  double peakval = hpost->GetBinCenter(MaxBin);

  // The total integral of the posterior
  double integral = hpost->Integral();

  // Keep count of how much area we're covering
  double sum = 0.0;

  // Counter for current bin. Going to get the error but counting outwards in each direction from the center til we reach fraction of events corresponding to 1 sigma
  int CurrBin = MaxBin;
  while (sum / integral < 0.6827 / 2.0 && CurrBin < hpost->GetNbinsX() + 1)
  {
    sum += hpost->GetBinContent(CurrBin);
    CurrBin++;
  }
  double sigma_p = fabs(hpost->GetBinCenter(MaxBin) - hpost->GetBinCenter(CurrBin));
  // Reset the sum
  sum = 0.0;

  // Reset the bin counter
  CurrBin = MaxBin;
  // Counter for current bin, now going counting bins down from HPD bin
  while (sum / integral < 0.6827 / 2.0 && CurrBin >= 0)
  {
    sum += hpost->GetBinContent(CurrBin);
    CurrBin--;
  }
  double sigma_m = fabs(hpost->GetBinCenter(CurrBin) - hpost->GetBinCenter(MaxBin));

  // Now do the double sided HPD
  sum = 0.0;
  int LowBin = MaxBin - 1;
  int HighBin = MaxBin + 1;
  double LowCon = 0.0;
  double HighCon = 0.0;
  while (sum / integral < 0.6827 && LowBin >= 0 && HighBin < hpost->GetNbinsX() + 1)
  {

    // Get the slice
    LowCon = hpost->GetBinContent(LowBin);
    HighCon = hpost->GetBinContent(HighBin);

    // If we're on the last slice and the lower contour is larger than the upper
    if ((sum + LowCon + HighCon) / integral > 0.6827 && LowCon > HighCon)
    {
      sum += LowCon;
      break;
      // If we're on the last slice and the upper contour is larger than the lower
    }
    else if ((sum + LowCon + HighCon) / integral > 0.6827 && HighCon >= LowCon)
    {
      sum += HighCon;
      break;
    }
    else
    {
      sum += LowCon + HighCon;
    }

    LowBin--;
    HighBin++;
  }

  double sigma_hpd = 0.0;
  if (LowCon > HighCon)
  {
    sigma_hpd = fabs(hpost->GetBinCenter(LowBin) - hpost->GetBinCenter(MaxBin));
  }
  else
  {
    sigma_hpd = fabs(hpost->GetBinCenter(HighBin) - hpost->GetBinCenter(MaxBin));
  }

  central = peakval;
  error = sigma_hpd;
  error_pos = sigma_p;
  error_neg = sigma_m;
}

// **************************
// Get the mean and RMS of a 1D posterior
void GetArithmetic(TH1D *const hpost, double &mean, double &error)
{
  // **************************
  mean = hpost->GetMean();
  error = hpost->GetRMS();
}

// **************************
// Get Gaussian characteristics
void GetGaussian(TH1D *&hpost, TF1 *&gauss, double &central, double &error)
{
  // **************************

  double mean = hpost->GetMean();
  double err = hpost->GetRMS();
  double peakval = hpost->GetBinCenter(hpost->GetMaximumBin());

  // Set the range for the Gaussian fit
  gauss->SetRange(mean - 1.5 * err, mean + 1.5 * err);
  // Set the starting parameters close to RMS and peaks of the histograms
  gauss->SetParameters(hpost->GetMaximum() * err * sqrt(2 * TMath::Pi()), peakval, err);

  // Perform the fit
  hpost->Fit(gauss->GetName(), "Rq");
  hpost->SetStats(0);

  central = gauss->GetParameter(1);
  error = gauss->GetParameter(2);
}

// **************************
// Function to get limits for all parameters from the input
void LoadInputVals()
{
  // **************************

  // Get asimov rates in a map
  TFile *asmvRatesFile = new TFile(asmvDistsFileName.c_str(), "READ");
  asmvRatesFile->cd();
  // Can't do GetObject for our own type (ParMap), so use a temporary std::map and cast
  std::map<std::string, double> *tempMap;
  asmvRatesFile->GetObject("AsimovRates", tempMap);
  asimovRates = (ParMap)*tempMap;

  // read fit config
  FitConfig mcConfig;
  FitConfigLoader mcLoader(mcmcConfigFile);
  mcConfig = mcLoader.LoadActive();

  constrMeans = mcConfig.GetConstrMeans();
  constrSigmas = mcConfig.GetConstrSigmas();
  mins = mcConfig.GetMinima();
  maxs = mcConfig.GetMaxima();

  ParMap tempMeans;
  ParMap tempSigmas;
  ParMap tempMins;
  ParMap tempMaxs;
  ParMap tempRates;

  // ttrees don't like hyphens in names so loop over parameters, switch any - to _ (currently only pmt-bg), and then set values for temporary ParMaps
  for (ParameterDict::iterator it = asimovRates.begin(); it != asimovRates.end(); ++it)
  {
    std::string tempName = it->first;
    if (tempName.find("-") != std::string::npos)
      tempName.replace(tempName.find("-"), 1, "_");
    tempMins[tempName] = mins[it->first];
    tempMaxs[tempName] = maxs[it->first];
    tempSigmas[tempName] = constrSigmas[it->first];
    tempMeans[tempName] = constrMeans[it->first];
    tempRates[tempName] = asimovRates[it->first];
  }

  // and now set all the maps to the temporary maps
  constrMeans = tempMeans;
  constrSigmas = tempSigmas;
  mins = tempMins;
  maxs = tempMaxs;
  asimovRates = tempRates;
}

// **************************
// Function to get limits for a single parameter from the already loaded input
void GetParLimits(std::string paramName, double &central, double &prior, double &down_error, double &up_error)
{
  // **************************

  central = 0.0; // Default to zero for anything not set in asimov rates file
  double error = 0.0;
  double sigmas = 1.0;

  if (asimovRates.find(paramName) != asimovRates.end())
  {
    central = asimovRates[paramName];
    prior = constrMeans[paramName];
    error = constrSigmas[paramName];
  }

  // If we have a prior, we want the uncertainty from the prior central value, not the asimov central value
  double x;
  if (prior)
    x = prior;
  else
    x = central;

  // We might be passed the valid range of parameter
  // Do a check to see if this is true
  if (x - error < mins[paramName])
  {
    down_error = mins[paramName];
  }
  else
  {
    down_error = x - error;
  }

  if (x + error > maxs[paramName])
  {
    up_error = maxs[paramName];
  }
  else
  {
    up_error = x + error;
  }

  if (prior)
    central = prior;
}

// *****************************
// Make the prefit plots
TH1D *MakePrefit(int nPar)
{
  // *****************************

  // Initialise a 1D plot over all parameters to 0
  TH1D *PreFitPlot = new TH1D("Prefit", "Prefit", nPar, 0, nPar);
  for (int i = 0; i < PreFitPlot->GetNbinsX() + 1; ++i)
  {
    PreFitPlot->SetBinContent(i + 1, 0);
    PreFitPlot->SetBinError(i + 1, 0.01); // If this is 0 root might not plot it
  }

  int count = 0;
  // You could have made asimov rates for more events than you fit over, so iterate over one of the maps from the fit config so you only get the fitted parameters
  // for(ParameterDict::iterator it = mins.begin(); it != mins.end(); ++it){
  for (int i = 0; i < bnames.size(); i++)
  {

    if (asimovRates[bnames.at(i).Data()] != 0)
    {
      PreFitPlot->SetBinContent(count + 1, asimovRates[bnames.at(i).Data()] / asimovRates[bnames.at(i).Data()]);
      PreFitPlot->GetXaxis()->SetBinLabel(count + 1, bnames.at(i));
      if (constrMeans[bnames.at(i).Data()] != 0)
      {
        PreFitPlot->SetBinContent(count + 1, constrMeans[bnames.at(i).Data()] / asimovRates[bnames.at(i).Data()]);
        PreFitPlot->SetBinError(count + 1, constrSigmas[bnames.at(i).Data()] / asimovRates[bnames.at(i).Data()]);
      }
    }
    // If the asimov rate is 0, don't normalise, but set to 1 so it's still represented as a ratio to the asimov rate
    else
    {
      PreFitPlot->SetBinContent(count + 1, asimovRates[bnames.at(i).Data()] + 1);
      PreFitPlot->GetXaxis()->SetBinLabel(count + 1, bnames.at(i));
      if (constrMeans[bnames.at(i).Data()] != 0)
      {
        PreFitPlot->SetBinContent(count + 1, constrMeans[bnames.at(i).Data()] + 1);
        PreFitPlot->SetBinError(count + 1, constrSigmas[bnames.at(i).Data()]);
      }
    }
    // Keep count of number of parameters to get the right bin each time. This depends on each iteration over a map to be done in the same order, which I think is ok
    count++;
  }

  PreFitPlot->GetYaxis()->SetTitle("Number of Events Relative to Asimov");

  PreFitPlot->SetDirectory(0);

  PreFitPlot->SetFillStyle(1001);
  PreFitPlot->SetFillColor(kBlue);
  PreFitPlot->SetMarkerStyle(21);
  PreFitPlot->SetMarkerSize(1.0);
  PreFitPlot->SetMarkerColor(kWhite);
  PreFitPlot->SetLineColor(kBlue);

  PreFitPlot->GetXaxis()->LabelsOption("v");

  return PreFitPlot;
}

// *****************************
// Set palette to be blue->red centered on white. Useful for the 2D ratio and difference plots
// if z axis is centered on 1 (for ratio) or 0 (for difference) as you can easily distinguish
// positive and negative changes
void BlueRedPalette()
{
  // *****************************
  // Take away the stat box
  gStyle->SetOptStat(0);
  // Make pretty correlation colors (red to blue)
  const int NRGBs = 5;
  TColor::InitializeColors();
  Double_t stops[NRGBs] = {0.00, 0.25, 0.50, 0.75, 1.00};
  Double_t red[NRGBs] = {0.00, 0.25, 1.00, 1.00, 0.50};
  Double_t green[NRGBs] = {0.00, 0.25, 1.00, 0.25, 0.00};
  Double_t blue[NRGBs] = {0.50, 1.00, 1.00, 0.25, 0.00};
  TColor::CreateGradientColorTable(5, stops, red, green, blue, 255);
  gStyle->SetNumberContours(255);
}

// *****************************
// Get the step number of the step with the maximum LLH
int GetMaxLLHStep(TTree *chain)
{
  // *****************************

  int maxLLHStep = 0;
  double LLH = 0;
  double maxLLh = -999;
  chain->SetBranchAddress("LogL", &LLH);

  // Loop over steps
  for (int i = 0; i < chain->GetEntries(); i++)
  {
    // Get step
    chain->GetEntry(i);

    // update max llh found so far, and the current highest step number
    if (LLH > maxLLh)
    {
      maxLLh = LLH;
      maxLLHStep = i;
    }
  }

  return maxLLHStep;
}

// *****************************
// Get the parameter value at a given step. Intention is for step to be the one with the
// maximum LLH in the chain, so you find the parameter value there. But could be used
// for any step (i.e maxStep could be any step number)
double GetMaxLLHPar(TTree *chain, int maxStep, std::string parname)
{
  // *****************************
  double parval = 0;

  chain->SetBranchAddress(parname.c_str(), &parval);
  chain->GetEntry(maxStep);
  return parval;
}

// *****************************
// Make a canvas comparing two 1D plots with ratio panel
TCanvas *CompareProjections(TH1D *prefit, TH1D *postfit)
{
  // *****************************

  prefit->Sumw2();
  postfit->Sumw2();

  postfit->GetXaxis()->SetTitleSize(0.07);
  postfit->GetXaxis()->SetLabelSize(0.06);
  postfit->GetYaxis()->SetTitleSize(0.07);
  postfit->GetYaxis()->SetLabelSize(0.06);

  // copy histograms for dividing to get ratio
  TH1D *postfitn = (TH1D *)postfit->Clone("postfitn");
  TH1D *prefitn = (TH1D *)prefit->Clone("prefitn");
  prefitn->Divide(prefit);
  postfitn->Divide(prefit);

  TCanvas *c2 = new TCanvas("c2", "c2", 0, 0, 1500, 1000);
  c2->SetGrid();

  postfitn->GetYaxis()->SetRangeUser(0.989, 1.011);
  postfitn->GetYaxis()->SetTitle("Ratio");
  postfitn->SetTitle("");

  // Make one canvas with both plots
  double uppermin = 0.35;
  double middlemin = 0.237;
  TPad *lower = new TPad("lower", "pad", 0, 0, 1, middlemin);
  lower->SetGrid();
  TPad *upper = new TPad("upper", "pad", 0, middlemin, 1, 1);
  upper->SetGrid();
  upper->SetBottomMargin(0.01);
  lower->SetTopMargin(0.01);
  lower->SetBottomMargin(0.58);
  upper->Draw();
  lower->Draw();
  c2->cd();

  // divided
  lower->cd();

  postfitn->SetTitle("");
  postfitn->GetYaxis()->SetTitle("Ratio");
  postfitn->GetYaxis()->CenterTitle();
  postfitn->GetYaxis()->SetTitleSize(0.16);
  postfitn->GetYaxis()->SetTitleOffset(0.32);
  postfitn->GetYaxis()->SetNdivisions(305, kTRUE);
  postfitn->GetYaxis()->SetLabelSize(0.1);
  postfitn->GetXaxis()->SetTitleSize(0.17);
  postfitn->GetXaxis()->SetTitleOffset(1.4);
  postfitn->GetXaxis()->SetLabelSize(0.13);
  postfitn->SetLineColor(kBlue);
  prefitn->SetLineColor(kRed);
  postfitn->Draw("L hist");
  prefitn->Draw("L hist same");

  // Draw line at unity for ratio plot
  TLine *line = new TLine(postfitn->GetXaxis()->GetBinLowEdge(1), 1.0, postfitn->GetXaxis()->GetBinUpEdge(postfitn->GetXaxis()->GetNbins()), 1.0);
  line->SetLineWidth(1.5);
  line->SetLineColor(kBlack);
  line->SetLineStyle(2);
  line->Draw();

  // un-divided, not- bin normalised
  upper->cd();
  postfit->SetTitle("");
  postfit->GetYaxis()->SetTitle("Events");
  postfit->GetYaxis()->SetTitleOffset(1.015);
  postfit->GetYaxis()->SetTitleSize(0.05);
  postfit->GetYaxis()->SetLabelSize(0.04);
  postfit->SetLineColor(kBlue);
  prefit->SetLineColor(kRed);
  postfit->Draw("e1");
  prefit->Draw("e1 same");

  // Draw legend
  TLegend *l = new TLegend(0.15, 0.7, 0.3, 0.85);
  l->AddEntry(prefit, "Prefit", "l");
  l->AddEntry(postfit, "Postfit", "l");
  l->SetFillColor(0);
  l->Draw();

  return c2;
}
