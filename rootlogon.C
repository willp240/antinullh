
// The rootlogon
void rootlogon()
{

  TStyle* snoStyle = new TStyle("snoplus", "SNO+ plots style for publications");

  // Use plain black on white colors
  snoStyle->SetFrameBorderMode(0);
  snoStyle->SetCanvasBorderMode(0);
  snoStyle->SetPadBorderMode(0);
  snoStyle->SetPadBorderSize(0);
  snoStyle->SetPadColor(0);
  snoStyle->SetCanvasColor(0);
  snoStyle->SetTitleColor(0);
  snoStyle->SetStatColor(0);
  // snoStyle->SetFillColor(0); // needs to be commented out so that TPalette works on TH2

  // Use bold lines
  snoStyle->SetHistLineWidth(2);
  snoStyle->SetLineWidth(2);

  // No title, stats box or fit as default
  snoStyle->SetOptTitle(0);
  snoStyle->SetOptStat(0);
  snoStyle->SetOptFit(0);

  // Postscript dashes
  snoStyle->SetLineStyleString(2, "[12 12]"); // postscript dashes

  // Text style and size
  snoStyle->SetLabelOffset(0.01, "x");
  snoStyle->SetTickLength(0.015, "x");
  snoStyle->SetTitleOffset(1.0, "x");
  snoStyle->SetLabelOffset(0.01, "y");
  snoStyle->SetTickLength(0.015, "y");
  snoStyle->SetTitleOffset(0.8, "y");
  snoStyle->SetLabelOffset(0.01, "z");
  snoStyle->SetTickLength(0.015, "z");
  snoStyle->SetTitleOffset(0.7, "z");
  snoStyle->SetLabelFont(132, "x");
  snoStyle->SetLabelFont(132, "y");
  snoStyle->SetLabelFont(132, "z");
  snoStyle->SetTitleFont(132, "x");
  snoStyle->SetTitleFont(132, "y");
  snoStyle->SetTitleFont(132, "z");
  snoStyle->SetLabelSize(0.05, "x");
  snoStyle->SetTitleSize(0.06, "x");
  snoStyle->SetTitleColor(1, "x");
  snoStyle->SetLabelSize(0.05, "y");
  snoStyle->SetTitleSize(0.06, "y");
  snoStyle->SetTitleColor(1, "y");
  snoStyle->SetLabelSize(0.05, "z");
  snoStyle->SetTitleSize(0.06, "z");
  snoStyle->SetTitleColor(1, "z");
  snoStyle->SetPadTickX(1);
  snoStyle->SetPadTickY(1);

  // Legends
  snoStyle->SetLegendBorderSize(0);
  snoStyle->SetLegendFont(132);
  snoStyle->SetLegendFillColor(0);

  // Graphs - set default marker to filled square
  snoStyle->SetMarkerStyle(21);

  // Palette for TH2 plots
  // Please try to use a Colour Vision Deficiency (CVD) friendly palettes.
  // You can make your own or use one of the ROOT pre-defined palettes on https://root.cern.ch/doc/master/classTColor.html
  snoStyle->SetPalette(kInvertedDarkBodyRadiator);

  // SNO+ Preliminary label
  snoStyle->SetTextFont(132);
  snoStyle->SetTextSize(0.06);
  snoStyle->SetTextAlign(32); // Right (horizontally) and center (vertically) adjusted

  gROOT->SetStyle("snoplus");
}