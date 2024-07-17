import ROOT
import glob
import os

# do these ones
smooth_names = ["tl208_id", "tl208", "bi214_id", "bi214_od"]
just_do = []

smooth_all_below = 2

def combine_projections(old, *hists):
    new = old.Clone()
    new.Scale(0)
    for ihist, hist in enumerate(hists):
        rbin = ihist + 1
        for ebin in range(1, hist.GetNbinsX()+ 1):
            new.SetBinContent(new.GetBin(ebin, rbin), hist.GetBinContent(ebin))
    return new

pdfs = glob.glob("../originals/*tl208*.root") + \
       glob.glob("../originals/*bi214*.root") + \
       glob.glob("../originals/*pmt*.root")

for x in twod_pdfs:
    can = ROOT.TCanvas()
    name = os.path.splitext(os.path.basename(x))[0]

    # want to do it to the twod_pdfs, not the projections
    if "proj" in name:
        continue
   
    f = ROOT.TFile(x)
    h = f.Get(f.GetListOfKeys()[0].GetName())

    if(len(just_do) != 0 and (name not in just_do)):
        print "skipping", name
        continue

    smoothing = (name in smooth_names) 

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
        projs.append(proj)

    combine_projections(h, *projs).SaveAs(os.path.basename(x))
