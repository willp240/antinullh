# Usage: python3 plotFixedOscLLHWhole -i <fitresulttree.root>
# Output OscLLH.png & OscLLH.pdf which record a 2D LLH with its profilellh distribution

import ROOT
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import os
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import matplotlib.font_manager as fm
import argparse
def SNOPLUSSTYLE():
    #personal_path = os.path.dirname(os.path.realpath(__file__))
    personal_path = "~/Library/Times_New_Roman_Normal.ttf"
    font_file = 'Times_New_Roman_Normal.ttf'
    font_path = personal_path + '/' + font_file
    plt.rcParams['axes.labelsize'] = 24
    # the size here only affects the tick labels because everything else is set for individual plots
    paper_font = fm.FontProperties(fname=font_path, size = 24)
    plt.rcParams['mathtext.default'] = 'regular'

    # set position of axis labels
    plt.rcParams["xaxis.labellocation"] = 'right'
    plt.rcParams["yaxis.labellocation"] = 'top'

    # set global parameters with rcParams -- copy whole block below before plotting code to set plotting template globally
    plt.rcParams['xtick.labelsize'] = 20   # font size for x-axis tick labels
    plt.rcParams['ytick.labelsize'] = 20
    # set height and width of big markings on axis x
    plt.rcParams['xtick.major.size'] = 8
    plt.rcParams['xtick.major.width'] = 1.6
    # set height and width of small markings on axis x
    plt.rcParams['xtick.minor.size'] = 3
    plt.rcParams['xtick.minor.width'] = 1.6
    # set height and width of big markings on axis y
    plt.rcParams['ytick.major.size'] = 8
    plt.rcParams['ytick.major.width'] = 1.6
    # set height and width of small markings on axis y
    plt.rcParams['ytick.minor.size'] = 3
    plt.rcParams['ytick.minor.width'] = 1.6
    # set thickness of axes
    plt.rcParams['axes.linewidth'] = 1.6
    # set plot background color
    plt.rcParams['figure.facecolor'] = 'white'
    # set plot aspect ratio -- change according to needs
    plt.rcParams['figure.figsize'] = (8.5, 6.5)
    # set padding (between ticks and axis label)
    plt.rcParams['xtick.major.pad'] = '6'
    plt.rcParams['ytick.major.pad'] = '6'
    # set padding (between plot and title)
    plt.rcParams['axes.titlepad'] = 12
    # set markings on axis to show on the inside of the plot, can change if needed
    plt.rcParams['xtick.direction'] = 'in'
    plt.rcParams['ytick.direction'] = 'in'

    # set ticks on both sides of x and y axes
    plt.rcParams['xtick.top'] = True
    plt.rcParams['xtick.bottom'] = True
    plt.rcParams['ytick.left'] = True
    plt.rcParams['ytick.right'] = True

    # define custom styles -- add **set_style to plot, as shown in EXAMPLE PLOT 1 below

    # each of these are pretty self-explanatory, can modify them if needed e.g. can set custom markers and colors if different data types are plotted

    histogram_style = {
        'histtype': 'step', 
        'color': 'blue',
        'alpha': 0.7,
        'linewidth': 2
    }

    scatter_style = {
        'marker': 's',
        'color': 'black',
        's': 25
    }

    errorbar_style = {
        'linestyle': 'None',
        'color': 'black',
        'capsize': 1.5
    }

    line_plot_style = {
        'linestyle': '-',
        'color': 'red',
        'linewidth': 2.5
    }
# Format the y-axis tick labels to be scaled
def scaled_formatter(x, pos):
    #The two args are the value and tick position
    #print('%1.0f' % (x/1e-5))
    return '%1.0f' % (x/1e-5)

if __name__ == "__main__":
    parser = argparse.ArgumentParser("Submit jobs to condor")
    parser.add_argument("-i", "--i", type=str, default="", help='fit_result_tree.root path')
    args = parser.parse_args()

    input_dir = os.path.dirname(args.i)
    output_dir = os.path.join(input_dir, "plots")
    os.makedirs(output_dir, exist_ok=True)
    # Load ROOT file and tree
    file = ROOT.TFile.Open(args.i)
    tree = file.Get("fitResults")
    nEntries = tree.GetEntries()

    # Extract tree data
    theta_list, deltam_list, llh_list = [], [], []
    for i in range(nEntries):
        tree.GetEntry(i)
        theta_list.append(tree.theta12)
        deltam_list.append(tree.deltam21)
        llh_list.append(tree.LLH)

    theta = np.array(theta_list)
    deltam = np.array(deltam_list)
    llh = np.array(llh_list)

    # Calculate 2ΔlnL
    llh_min = np.min(llh)
    llh_2delta = 2 * (llh - llh_min)

    # Binning for 2D histogram
    nbins = int(np.sqrt(nEntries))
    theta_bins = np.linspace(theta.min(), theta.max(), nbins + 1)
    deltam_bins = np.linspace(deltam.min(), deltam.max(), nbins + 1)

    # Histogram 2D
    hist2d, xedges, yedges = np.histogram2d(
        theta, deltam, bins=[theta_bins, deltam_bins], weights=llh_2delta
    )
    counts, _, _ = np.histogram2d(theta, deltam, bins=[theta_bins, deltam_bins])
    # element wise: gist2d/counts, and avoiding zero in counts(If counts = 0, output = 0)
    hist2d_avg = np.divide(hist2d, counts, out=np.zeros_like(hist2d), where=counts != 0)

    # Profiles
    profile_theta = np.min(hist2d_avg, axis=1)  # Min over deltam
    profile_deltam = np.min(hist2d_avg, axis=0)   # Min over theta
    theta_centers = 0.5 * (xedges[:-1] + xedges[1:])
    deltam_centers = 0.5 * (yedges[:-1] + yedges[1:])
    #print(profile_theta-profile_deltam)

    # Plot setup
    SNOPLUSSTYLE()
    fig = plt.figure(figsize=(10, 8))
    gs = gridspec.GridSpec(2, 2, width_ratios=[3, 1], height_ratios=[1, 3], hspace=0.0, wspace=0.0)

    # Share x and y axes as needed
    ax_heat = fig.add_subplot(gs[1, 0]) #2D LLH
    ax_top = fig.add_subplot(gs[0, 0]) # theta12 
    ax_right = fig.add_subplot(gs[1, 1]) # delm12
    cbar_ax = fig.add_subplot(gs[0, 1])  #colorbar

    # Setting Minor Tick On
    ax_heat.minorticks_on()
    ax_top.minorticks_on()
    ax_right.minorticks_on()
    cbar_ax.minorticks_on()
    # Top left: θ12 profile
    ax_top.plot(theta_centers, profile_theta, color='black')
    ax_top.set_xlim(theta_centers.min(), theta_centers.max())
    ax_right.set_label(r"$2\,\Delta\ln\mathcal{L}$")
    ax_top.set_xticklabels([])

    # Bottom left: 2D LLH plot
    X, Y = np.meshgrid(xedges, yedges)
    ax_heat.yaxis.set_major_formatter(mticker.FuncFormatter(scaled_formatter))
    ax_heat.set_ylabel(r"$\Delta m_{21}^2$ ($\times 10^{-5}$ MeV$^2$)")
    pc = ax_heat.pcolormesh(X, Y, hist2d_avg.T, shading='auto', cmap='viridis')
    ax_heat.set_xlabel(r"$\theta_{12}$")
    # plot contour for 1 and 2 sigma
    # Contour levels: 1σ and 2σ
    contour_levels = [2.30, 6.18] 
    CS = ax_heat.contour(
        theta_centers, deltam_centers, hist2d_avg.T,  # Note: need centers, not bin edges
        levels=contour_levels,
        colors=['white', 'white'],
        linewidths=2,
        linestyles=['solid', 'dashed']
    )
    fmt = {2.30: r"$1\sigma$", 6.18: r"$2\sigma$"}
    ax_heat.clabel(CS, fmt=fmt, inline=True, colors='black')

    # Bottom right: deltam profile
    ax_right.plot(profile_deltam, deltam_centers, color='black')
    ax_right.set_ylim(deltam_centers.min(), deltam_centers.max())
    ax_right.set_label(r"$2\,\Delta\ln\mathcal{L}$")
    ax_right.set_yticklabels([])

    # Top right: colorbar
    fig.colorbar(pc, cax=cbar_ax, orientation='vertical')
    cbar_ax.yaxis.set_ticks_position('right')
    cbar_ax.set_label(r"$2\,\Delta\ln\mathcal{L}$")

    #plt.tight_layout()
    #plt.show()
    plt.savefig(os.path.join(output_dir, "OscLLH.png"))
    plt.savefig(os.path.join(output_dir, "OscLLH.pdf"))
    print(f"successfully saving .png and .pdf to {output_dir}")
