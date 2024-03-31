import ROOT
import glob
import os

COLORS = [ROOT.kRed-4, ROOT.kBlue-4, ROOT.kGreen+2, ROOT.kOrange+1,
          ROOT.kMagenta+1, ROOT.kCyan+1, ROOT.kOrange-2, ROOT.kRed+2,
          ROOT.kBlue-9, ROOT.kSpring+4, ROOT.kGray+1, ROOT.kPink+1,
          ROOT.kViolet+6]

# do these ones
smooth_names = ["tl208_id", "tl208", "bi214_id", "bi214_od"]
just_do = []


smooth_all_below = 2


def combine_projections(old, *hists):
    new = old.Clone()
    new.Scale(0)
    for ihist, hist in enumerate(hists):
        rbin = ihist + 1
        for ebin in range(1, hist.GetNbinsX() + 1):
            new.SetBinContent(new.GetBin(ebin, rbin), hist.GetBinContent(ebin))
    return new

twod_pdfs = glob.glob("../originals/*tl208*.root") + \
            glob.glob("../originals/*bi214*.root") + \
            glob.glob("../originals/*pmt*.root")

for x in twod_pdfs:
    can = ROOT.TCanvas()
    name = os.path.splitext(os.path.basename(x))[0]

    # want to do it to the twod_pdfs, not the projections
    if "proj" in name:
        continue

    if name == "tl208" or name == "tl208_id":
        leg = ROOT.TLegend(0.15, 0.5, 0.5, 0.9)
    else:
        leg = ROOT.TLegend(0.55, 0.5, 0.9, 0.9)
   
    f = ROOT.TFile(x)
    h = f.Get(f.GetListOfKeys()[0].GetName())

    if(len(just_do) != 0 and (name not in just_do)):
        print "skipping", name
        continue

    smoothing = (name in smooth_names) 
    if smoothing:
        name += "_smoothed"

    projs = []
    for i in xrange(1, h.GetNbinsY() + 1):
        # make a copy of the projection
        proj = h.ProjectionX("bin{0}".format(i),i, i).Clone()
        before = proj.Integral()

        if smoothing or i <= smooth_all_below:
            proj.Smooth(1)

        # make sure the integral is preserved
        after = proj.Integral()
        if after != 0:
            proj.Scale(before/after)

        proj.SetStats(0)
        proj.Draw("SAME")

        can.Update()
        leg.AddEntry(proj, "bin{0}".format(i-1))
        proj.GetXaxis().SetTitle("Fit Energy/MeV")
        proj.GetYaxis().SetTitle("Probability")
        proj.SetTitle("")
        proj.SetLineColor(COLORS[i])
        projs.append(proj)
    leg.Draw("same")
    combine_projections(h, *projs).SaveAs(os.path.basename(x))
    print os.path.basename(x)
    can.SaveAs("{0}.pdf".format(name))
