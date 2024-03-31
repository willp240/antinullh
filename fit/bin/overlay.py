import ROOT
import ConfigParser

class Grouper:
    def __init__(self):
        self.hists = {}

    def add(self, hist, group):
        if hist.Integral() == 0 or hist.Integral() != hist.Integral(): 
            return
        try:
            self.hists[group].Add(hist)
        except KeyError:
            self.hists[group] = hist

def assign_group(parser, name):
    try:
        return parser.get(name, "plot_group")
    except ConfigParser.NoOptionError as e:
        return name

def grab_hist(filename):
    f = ROOT.TFile(filename)
    h = f.Get(f.GetListOfKeys()[0].GetName())
    h.SetDirectory(0)    
    return h

def stack(order, hists, labels, colors):
    stack = ROOT.THStack("fit", "")
    leg   = ROOT.TLegend(0.75, 0.65, 0.95, 0.95)

    for name in order:
        try:
            stack.Add(hists[name])
        except KeyError as e:
            print "Warning: found no hists for group " + name
            continue
        t_col = ROOT.TColor.GetColor("#{0}".format(colors[name]))
        hists[name].SetLineColor(t_col)
        hists[name].SetFillColor(t_col)
        
    for name in reversed(order):
        try:
            leg.AddEntry(hists[name], labels[name], "FL")
        except KeyError as e:
            pass
    return stack, leg

def stack_slice(order, hists, labels, colors, bin):
    stack = ROOT.THStack("fit_{0}".format(bin), "")
    leg   = ROOT.TLegend(0.65, 0.55, 0.88, 0.88)
    projs = {}
    for name in order:
        try:  
            proj = hists[name].ProjectionX("{0}_{1}slice_px".format(name, bin), bin, bin)
            proj.SetDirectory(0)
            projs[name] = proj
            stack.Add(proj)
        except KeyError as e:
            print "Warning: found no hists for group " + name
            continue
        t_col = ROOT.TColor.GetColor("#{0}".format(colors[name]))
        projs[name].SetLineColor(t_col)
        projs[name].SetFillColor(t_col)

    for name in reversed(order):
        try:
            if (name=="alphan"):
                leg.AddEntry(projs[name], "(#alpha, n)", "FL")
                continue
            if (name=="2v"):
                leg.AddEntry(projs[name], "2#nu#beta#beta", "FL")
                continue
            if (name=="0v"):
                leg.AddEntry(projs[name], "0#nu#beta#beta", "FL")
                continue
            if (name=="b8_nue"):
                leg.AddEntry(projs[name], "^{8}B #nu", "FL")
                continue
            if (name=="U Chain"):
                leg.AddEntry(projs[name], "^{238}U Chain", "FL")
                continue
            if (name=="Th Chain"):
                leg.AddEntry(projs[name], "^{232}Th Chain", "FL")
                continue
            leg.AddEntry(projs[name], labels[name], "FL")            
        except KeyError as e:
            pass
    return stack, leg, projs


def savestack(data, stack_, leg, name, x_title, y_title, title, result_dir):

   #ROOT.gROOT.SetStyle("thesisStyleMinipage")
   #ROOT.gROOT.ForceStyle()
    
    
    # draw it
    can =  ROOT.TCanvas()
    ROOT.gPad.SetLogy()
    stack_.Draw()
    stack_.GetXaxis().SetTitle(x_title)
    stack_.GetYaxis().SetTitle(y_title)
    stack_.GetYaxis().SetRangeUser(0.5, 10000)

    stack_.GetYaxis().SetLimits(0.5, 10000)
    stack_.SetTitle(title)
    stack_.Draw()
    stack_.SetMinimum(0.5)
    stack_.SetMaximum(10000)
    can.Modified()

    data.SetMarkerSize(0)
    
    data.Draw("same PE")
    data.SaveAs(os.path.join(result_dir, name + "data.root"))
    leg.AddEntry(data, "Fake Data", "PE")
    leg.Draw("same")
    can.SaveAs(os.path.join(result_dir, name + "stacked_fit.root"))
    can.SaveAs(os.path.join(result_dir, name + "stacked_fit.pdf"))
    
    
    # now work out a chi square
    s_hist = stack_.GetStack().Last()
    chi_square = 0
    for i in xrange(1, s_hist.GetNbinsX() + 1):
        try:
            chi_square += (s_hist.GetBinContent(i) -  data.GetBinContent(i)) ** 2 / s_hist.GetBinContent(i)
        except ZeroDivisionError as e:
            print "zero probability bin! #", i, "  @ ", s_hist.GetXaxis().GetBinCenter(i)
    
            
    print "\n\nchisq/bin = ", chi_square, " / ", s_hist.GetNbinsX()
    print "p-value = ", ROOT.TMath.Prob(chi_square, s_hist.GetNbinsX())


    print "Integrals:\n"
    print "\t Fake Data : ", data.Integral()
    print "\t Stack : ", s_hist.Integral()
    print "\t Sigma  : ", (data.Integral() - s_hist.Integral())/math.sqrt(s_hist.Integral())

    # now calculate the diff
    s_hist.Sumw2()
    data.Add(s_hist, -1)
    data.SaveAs(os.path.join(result_dir, name + "fit_diff.root"))


if __name__ == "__main__":
    import argparse
    import glob
    import os
    import math
    parser = argparse.ArgumentParser()
    parser.add_argument("plot_config", type = str)
    parser.add_argument("event_config", type = str)
    parser.add_argument("result_dir", type=str)
    parser.add_argument("--proj_titles_file", type=str)
    args = parser.parse_args()
    data  = grab_hist(os.path.join(args.result_dir, "scaled_dists", "data.root"))
    data.Sumw2()

    # read 
    scaled_dists_p = [ x for x in glob.glob(os.path.join(args.result_dir, "scaled_dists", "*.root")) if not x.endswith("data.root")]
    names = [os.path.basename(x).split(".root")[0] for x in scaled_dists_p]
    
    
    scaled_dists = {}
    for name, path in zip(names, scaled_dists_p):
        scaled_dists[name] = grab_hist(path)

    # read from the config file if these guys are in a group
    parser = ConfigParser.ConfigParser()
    parser.read(args.event_config)
    groups = dict((x, assign_group(parser, x)) for x in names)
    tex_labels = dict((x, parser.get(x, "tex_label")) for x in names)

    # filter them into their groups
    grouper = Grouper()
    for name in names:
        p = scaled_dists[name].ProjectionY("blah", 1, 1)
        grouper.add(scaled_dists[name], groups[name])

    # read the plot config
    parser = ConfigParser.ConfigParser()
    parser.read(args.plot_config)

    for gp in grouper.hists.keys():
        try:
            tex_labels[gp] = parser.get(gp, "tex_label")
        except ConfigParser.NoOptionError as e:
            if not gp in tex_labels.keys():
                tex_labels[gp] = gp

        except ConfigParser.NoSectionError as e:
            pass
    colors = {}
    for name in grouper.hists.keys():
        try:
            colors[name] =  parser.get(name, "color")
        except ConfigParser.NoSectionError as e:
            colors[name] = "aaaaaa"
        except ConfigParser.NoSectionError as e:
            colors[name] = "aaaaaa"

    order  = parser.get("summary", "plot_order").split(",")
    x_title  = parser.get("titles", "x_axis")
    y_title  = parser.get("titles", "y_axis")

    if all(type(x).__name__ == "TH1D" for x in scaled_dists.values()):
        stack_, leg = stack(order, grouper.hists, tex_labels, colors)
        savestack(data, stack_, leg, "", x_title, y_title, "", args.result_dir)

    elif all(type(x).__name__ == "TH2D" for x in scaled_dists.values()):
        n_y_bins = grouper.hists.values()[0].GetNbinsY()
        if args.proj_titles_file is not None:
            with open(args.proj_titles_file) as f:
                proj_titles = f.read().splitlines()

            if len(proj_titles) != n_y_bins:
                raise ValueError("If you specify the projection titles there should be one per y bin ({0}), you gave {1}".format(n_y_bins, len(proj_titles)))
        else:
            proj_titles = ["projection_{0}".format(iBin) for iBin in xrange(n_y_bins + 1)]

        # want to stack and save each slice independently
        for iBin in xrange(1, n_y_bins + 1):
            stack_, leg, projs = stack_slice(order, grouper.hists, tex_labels, 
                                             colors, iBin)            
            savestack(data.ProjectionX("_px", iBin, iBin, "e"), stack_, leg, "proj{0}_".format(iBin), x_title, y_title, proj_titles[iBin - 1], args.result_dir)
    else:
        print "Can't make sense of the histograms.. Only TH1D and TH2D are supported. not doing anything"
