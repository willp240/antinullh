<h1>Antinu Likelihood Analysis How To Guide</h1>

This file is a guide on how to run the most commonly used apps in the antinu oxo likelihood analysis.  

The first thing you'll want to do is get the environment variables in `env.sh` pointing to your OXO install, this `antinullh` directory, and the top level where data will be written to. From now on, paths in this document are relative to `antinullh` (except for postfit outputted files).

Most of the apps and source code live within the top directory. Within `./src`, there is code for interfacing with config files and some other useful classes. The config files themselves live in `./cfg`.

<h2>Config Loaders</h2> 
For each type of config file, there is a class defined in `./src`, along with a loader for that class. Each of these follows a similar structure: the config class has an attribute for each field in the config file, and the config loader classes read values in from the config files and return an instance of the corresponding config class. 

There are also other useful classes and files in `./src`:

<h2>Systematic Factory</h2> 
The systematic factory, when given a systematic name, type, and list of associated fit parameters, will return a systematic object to be used in the analysis. For some systematics, a function, oscillation grid, and reactor-distance map can also be handed to the systematic factory. The allowed types differ from just the OXO inherited systematic classes as they can also define the function used, as these are defined inline with lambdas to pass the oscillation grid and map if needed. Currently the allowed `types` are:
-`scale`: A linear Scale systematic
-`shift`: A Shift systematic
-`sqroot_scale_conv`: A Gaussian Convolution where the width of the Gaussian runs with the square root of the axis in question
-`scale_function`: A ScaleFunction systematic where the function used modifies Birk's Constant
-`shape`: A Shape systematic, where the function used modifies the oscillation probability

<h2>Functions</h2>
Functions that can be passed to systematics are defined in here. Currently they are Birk's Law and the Oscillation Probability calculations. However, to pass lambdas to these functions, they had to be defined inline so this file may soon be removed.

<h2>DistBuilder</h2>
The DistBuilder class can build a `BinnedED` from a combination of a PDF Config and a dataset. 

<h2>OscGrids</h2>
The oscillation probability calculation can be quite slow to do on the fly, so instead we can produce an "Oscillation Grid" (`OscGrid`). The probability of oscillation depends on the true neutrino energy, the distance travelled, and the oscillation parameters $\Delta m^2_{21}$ and $\theta_{12}$, so we can pre-calculate the probabilities for different values of each of those variables, and the linearly interpolate between them. As the distances to reactor cores are all fixed distinct values, rather than interpolate, we can produce a 3D grid for each reactor core distance.

The `OscGrid` class has attributes for the min and max values and number of grid points for each dimension (true neutrino energy, $\Delta m^2_{21}$ and $\theta_{12}$), along with a vector of probabilities where each element represents the probability at a single point on the grid. An `OscGrid` can be constructed using just the axis ranges and number of points, along with the filename the grid will be written to and the distance it applies for. You can then calculate the probabilities at each grid point using the `CalcGrid` method. This can be written to the output file with the `Write` method. Alternatively you can load up a pre-written and calculated grid with the `Load` method. You can then obtain the probability for a given set of parameters with the `Evaluate` method. There are also 'Getters' for the probability vector, and vectors containing the grid points for each dimension. The calculation of the probability uses code copied from the RAT-tools/AntinuTools repo. We will look at making that code portable so we can call it rather than copy it in the future. The interpolation for calculating probabilities between grid points is done using functions which are stolen from the OXO Solar Analysis (https://github.com/dcookman/solar_analysis/). 

<h2>Utilities</h2>
This file contains utility functions, mostly related to the Oscillation Grids. There is a function to turn the `OscGrids` max, min, and number of point values for a given dimension into a vector, along with a function to find the grid points that surround a given point in the `OscGrid`. These are both stolen from https://github.com/dcookman/solar_analysis. There is also a function here to load the reactor core indices and distances from SNO+ into a map from a JSON file.

<h2>Configs</h2>
There are several config files you'll need for running different apps:

<h3>Event</h3>
This file should contain information on the types of events you want to include in the analysis. An example is shown in `cfg/event_config.ini`. There should be a section, `summary`, which controls which event types are active and where the MC files live. Then there's a table for each event type, which contains the expected rate, a Latex name for axis titles, the exact location of files, the dimensions of the PDF, and the groups the PDF should belong to (different syestematics apply to different groups, see the OXO documentation for more info).

<h4>Summary</h4>
- `active`: The event types to be used in the analysis. These should be comma separated, and each should also be the name of an Event field also defined in the file. You can also use `all` to use all events defined in the file
- `inactive`: Event types defined in the config that aren't being used in the analysis at this time
- `orig_base_dir`: The directory where the original (unpruned) MC files have been downloaded to
- `pruned_ntup_dir`: The directory where the pruned MC files should be written to, and later read from
- `pdf_dir`: Where files containing histograms for each PDF should be written to

<h4>Event Type Tables</h4>
The name of these tables should be the name of the event type.
- `rate`: The expected rate of the events per year
- `tex_label`: A latex name for the event type
- `ntup_files`: File path for the original unpruned files, relative to `orig_base_dir`
- `dimensions`: Number of dimensions the PDF should have
- `groups`: The groups the PDFs are part of. Each systematic is applied to one group (or all event types)

<h3>Fit</h3>
This file should contain information on the parameters of the fit, and the fit itself. An example is shown in `cfg/fit_config.ini`.

<h4>Summary</h4>
-`datafile`: The file containing the data events
-`asimov`: Bool to determine if an Asimov or data fit should be run. If 1/True, the scaled PDFs are used as the data to run an Asimov fit
-`livetime`: The livetime in years, used to scale the PDFs correctly to match the data
-`iterations`: The number of steps in the Markov Chain
-`burn_in`: The number of steps rejected at the start of the Markov Chain. This can be changed post analysis for most purposes
-`fit_dists`: Which PDFs should be fit. Can use `all`
-`output_directory`: Where to write output files To
-`n_steps`: Number of steps per iteration for Hamiltonian MCMC
-`epsilon`: Hamiltonian step size
-`hmc_iterations`: Number of steps in Hamiltonian MCMC
-`hmc_burn_in`: The number of steps rejected at the start of the Hamiltonian Markov Chain. This can be changed post analysis for most purposes
-`sigma_scale`: Global scaling factor each step size gets multiplied by
-`beeston_barlow`: Bool to determine if we should use the Beeston-Barlow method to account for MC stats uncertainty. If 1/True, it is used

<h4>Fit Parameter Tables</h4>
The name of these tables should be the name of the fit parameter.
-`min`: The minimum value the parameter is allowed to take
-`max`: The maximum value the parameter is allowed to take
-`sig`: The relative width of the Gaussian used to propose new step values (each gets scaled by summary:sigma_scale). Often referred to as the step size
-`nbins`: Number of bins used for plotting the parameter
-`constraint_mean`: Mean of the prior constraint on the parameter
-`constraint_sig`: Width of the prior constraint on the parameter
-`nom`: Nominal value of the parameter. Not used for the PDF normalisation parameters, as these just have a nominal value of rate(post cuts)*livetime

A note on HMCMC: There is currently the functionality to run some MCMC, followed by some HMCMC starting from the best fit points of the MCMC. However, this didn't give any improvement in results or efficiency over the straight MCMC. The functionality is preserved in case we ever want to use it in the future, but generally we'll just run MCMC, with a notional 1 step HMCMC ran to avoid problems with files not being created when they are expected to.
It's likely in the future we'll just remove the functionality.

<h3>PDF</h3>
This file should contain information on the pdf axes. An example is shown in `cfg/pdf_config.ini`.

<h4>Summary<h4>
-`build_order`: The axes used for the PDFs. These will be defined in the same file as axis tables. The order is important: the dimensions field in `events` config sets the number of dimensions, N, for a given PDF. The axes used will be the first N axes here
-`data_axes`: The axes the data will have
-`pdf_dir`: Where the PDFs should be saved to as histograms

<h4>PDF Axis Tables<h4>
The name of these tables should be the name of the axis.
-`nbins`: The number of bins in the axis
-`min`: The minimum value on the axis
-`max`: The maximum value on the axis
-`branch_name`: The name of the branch in the ntuple
-`tex_name`: A latex name for the axis to be used in plot labels

<h3>Syst</h3>
This file should contain information on the systematics used in the analysis. An example is shown in `cfg/syst_config.ini`.

<h4>Summary<h4>
-`active`: Which systematics (defined in the systematics tables) will be used in the analysis

<h4>Systematic Tables<h4>
The name of these tables should be the name of the systematic.
-`param_names`: Comma separated list of the names of the fit parameters that control the systematic. These should all match a fit parameter defined in the fit config. If this is empty, the systematic name is assumed to be the name of the fit parameter
-`dist_obs`: These are the observables that the PDFs the systematic will apply to have
-`trans_obs`: These are the observables that the systematic will use to make the transformation
-`type`: The type of systematic being applied. This is not just the OXO derived systematic class, as here we also determine which function we will use for function-based systematics. The options are defined in the SystFactory (this will be made clearer in the future)
-`function`: The name of the function used if the systematic needs one. There is some redundancy here with type at the moment but we'll improve that soon
-`group`: This is the group of pdfs the systematic will be applied to. If this is left blank it will apply to all PDFs

<h3>Oscillation Grid</h3>
This file should contain information about the oscillation grids you want to use in the analysis, both for producing them, and reading them.

<h4>Summary<h4>
-`filename`: Where the oscillation grids should be written to and read from
-`reactorsjson`: Where the reactor distances should be read from
-`mine`: The minimum true neutrino energy of the grid
-`maxe`: The maximum true neutrino energy of the grid
-`numvalse`: The number of grid points on the energy axis
-`mindm21sq`: The minimum $\Delta m^{2}_{21}$ value of the grid
-`maxdm21sq`: The maximum $\Delta m^{2}_{21}$ value of the grid
-`numvalsdm21sq`: The number of grid points on the $\Delta m^{2}_{21}$  axis
-`minssqth12`: The minimum sin$^2 \theta_{12}$ value of the grid
-`maxssqth12`: The maximum sin$^2 \theta_{12}$value of the grid
-`numvalsssqth12`: The number of grid points on the sin$^2 \theta_{12}$ axis

These classes and config files, along with all the OXO classes, are brought together in various apps in the `antinullh` top directory. These are described below in the order you will probably want to use the:

<h2>Calculating Reactor Distances</h2>

To calculate the oscillation probability for an MC reactor IBD event, we need to know how far it has travelled, so we need to know the distance from the reactor it was produced in to SNO+. Rather than we do this calculation on the fly at every event for every iteration of a fit, we precalculate the distance for each reactor core. You can do this by running:

>> ./bin/make_reactor_json

This will loop over cores defined in the REACTORS RATDB table, and calculate the distance from it to SNO+. This gets saved in a JSON file, along with the reactor's name, and a unique integer assigned to the core. The integer is used as a numerical dimension for the reactor neutrino PDFs.

<h2>Making Oscillation Grids</h2>

To produce a single reactor core `OscGrid`, run:

>> ./bin/make_osc_grids oscgrid_config index_to_write

In the submitting batch jobs section we will discuss running over all reactor cores at once.

<h2>Pruning Trees</h2>

Raw SNO+ ntuples, although much more lightweight than RATDS files, still contain much more information than we need at analysis level. So instead of carrying around that deadweight for the whole analysis, we first take the time to prune out the branches we don't need and save what we do need in new files to be used.  

The `make_trees` app does this for us. Having set your filepaths in your `events` config file, simply run:  

>> ./bin/make_trees cfg/event_config.ini

This looks for ntuple files for each of the event types specified in the `event` config (you can select/deselect event types with the active/inactive fields) and produces a new set of files containing only the relevant trees for the analysis (energy, position, fit validity, pulse shape discrimination values, and true neutrino energy and reactor name for reactor neutrinos). Be warned, if using all event types, this will take some time! (of order several hours).

<h2>LLH Scans</h2>

Right, we've now pruned our trees and used them to build PDFs and a simulated dataset(s). In theory, we're now ready to launch a fit! But let's be steady Eddies, Cautious Carols, and Nervous Nigels, and make sure everything is behaving as intended. Running fits can be computationally and storage expensive, so it would be a shame to find out after running them that something had previously gone wrong.  

One of way doing this is with likelihood scans. The `llh_scan` app first builds the same test statistic that the fit will use, and we will hand it the same dataset we intend to hand the fitter. For a single floated parameter, it will scan over a range of values centred at the nominal value. For each of 150 steps, it will rebuild the dataset by scaling the PDFs and applying any systematics all set at their nominal values, except the parameter in question which is set to the value we've reached in the scan. This is compared to the actual dataset using the test statistic (probably the binned LLH). The LLH is saved at each step, and once we reach the end of the scan for a parameter, it is set back to its nominal value, and we repeat for the next parameter.  

For the 'true' Asimov dataset, these scans should all minimise at 1 (where the x axis is parameter value / nominal value) i.e., by changing a parameter, you can't produce a dataset more similar to the target dataset than by using the exact values used to produce the target dataset. For the other forms of Asimov dataset, most parameters will minimise close to 1, but probably not exactly. The event types with higher rates will be closer to 1 as changes in these will have a greater impact on the LLH, and the small fluctuations causing the differences between the PDFs and Asimov dataset will be smaller proportionally.  

To run the LLH scan:  

>> ./bin/llh_scan cfg/fit_config.ini cfg/event_config.ini cfg/pdf_config.ini cfg/syst_config.ini cfg/oscgrid.ini

A root file `llh_scan.root` will be saved in the output directory specified in the fit config. In the file will be a plot of LLH vs parameter value for each parameter.   

<h2>Running a Fit</h2>

Now the time has come to run a fit. First you want to make sure everything in you fit config is looking sensible. The min and max values, and sigmas for each event type are fairly well tuned, so it's advised you leave these alone unless you know what you're doing. Certainly, for your first fit testing out the code I would leave it as is. If you have good reason to change them, you probably don't need to be reading this guide! You can run a fit by:

>> ./bin/fit_dataset cfg/fit_config.ini cfg/event_config.ini cfg/pdf_config.ini cfg/syst_config.ini cfg/oscgrid.ini

<h2>Submitting Batch Jobs</h2>

For running full fits, it’s advisable to run multiple fits at once in parallel as you probably want around 1 million steps. Every batch system will be different, but for submitting to a Condor based queue system there is a script, `util/submitCondor.py`, based on one for submitting RAT jobs originally from Josie (I think).  

You can submit N jobs with:  

>> python utils/submitCondor.py exec_name output_dir -r /path/to/this/repo/ -e /path/to/env/file/ -f cfg/fit_config.ini -i cfg/event_config.ini  -p cfg/pdf_config.ini -s cfg/syst_config.ini -o cfg/oscgrid.ini -n N -w walltime

The environment file you supply here should be whatever you use to set up ROOT, python, GSL, etc., not `antinullh/env.sh` (which should be sourced before running the python submission). The output directory will be created, with several subdirectories:

`error`, `output` and `log` contain the console outputs and batch log. `sh` and `submit` contain the files used to submit these jobs.

You can run this script with any of the apps. If you're not running a fit, you probably don't need multiple simultaneous jobs so can just run with N=1. If you're running `make_osc_grid` with the submission script, it will automatically loop over every reactor core in the reactors JSON file and run `make_osc_grid` for each in parallel.

<h2> Postfit Analysis</h2>

NOTE: This section is a bit out of date because the code hasn't been brought in line with the rest of the antinullh repo yet. It will be updated when the code is!

Now you’re fit has run, the real fun starts! For each individual chain, you’ll have a number of files and subdirectories. In `1dlhproj` and `2dlhproj`, you have projections of the LLH for each parameter, and each combination of two parameters. In each case all other parameters are marginalised over. The burn-in steps, as set in the fit config, are automatically not included.  

`config_log.txt` and `config_log_hmc.txt` save the configs used (these are almost certainly the same as each other). These are also saved in the `cfg` directory if running the python submission script in `util` (see below). `auto_correlations.txt` calculates how correlated the LLH is to the LLH at a step a different number of steps previous. This should hopefully be close to 0 after a few hundred. `fit_results.txt` contains each parameter value for the maximum LLH step. `scaled_dists` contains the PDFs for each event type, scaled by the event type parameter value at the maximum LLH step in this chain. Also saved is the sum of these, with each systematic’s value at the maximum LLH step applied.  

There will also be a root file (`fit_name_i.root`) which contains a tree. Each entry in the tree represents a step in the Markov Chain, and the leaves are the values of each parameter, as well as the LLH, step time, and acceptance rate. This is more for MCMC chain diagnostics and debugging than obtaining physics results but is useful for checking the fit has worked as expected.  

To combine multiple parallel chains, `hadd` these root files. When plotting from the `hadded` file, you can cut off the burn-in using the `Step` branch (as entry will now run from 0 to total number of combined steps, not 0 to number of steps in each file).  

There’s an app for making useful plots from the combined tree. This can be run as:  

>> ./bin/make_plots output_tree.root asimov_dataset.root scaled_postfit_dist.root fit_config_file.cfg correlations traces

`output_tree.root` is either a single outputted file or the `hadded` combination. `asimov_dataset.root` is the ‘true’ Asimov dataset to get the Asimov rates to compare to. `scaled_postfit_dist.root` is the sum of the PDFs scaled by the maximum LLH step values. If you’ve `hadded` chains, you can find the chain with the highest LLH step (so you know which `scaled_dist` to use) using `util/findMaxStep.C`:  

>> root -l -b -q 'findMaxStep.C("dataset_directory_name", "fit_name", number_of_chains )'

`correlations` and `traces` determine if you draw the correlations between all combinations of two parameters (parameter 1 vs parameter 2 for all steps) and the trace of each parameter (parameter vs Step). These are both very time consuming. Leaving them empty means neither will run (note that the code just looks for a non-NULL string, so setting them to ‘0’ means they will be run!).  

The `make_plots` code is pretty messy, with features having been tacked on as time has gone on. But it does a job. The outputs are a root file of histograms, and a PDF with a bunch of histograms drawn to separate pages. For each parameter, you’ll get a LLH distribution marginalised over all other parameters. This is essentially a frequency histogram of parameter value for all steps. The central postfit value is calculated from this histogram in four different ways: the arithmetic mean of the histogram, the mean of a Gaussian fitted to it, the highest posterior density (mode), and the value at the maximum LLH step. These central values (and associated uncertainties) are combined into a single plot after the individual parameter plots. Prefit and postfit (using the maximum LLH step values) distributions are plotted and compared in one and two dimensions. If `traces` was set, a plot of parameter value vs `Step` is also drawn. If `correlations` are a run, a `colz` plot of parameter values at all steps are drawn.  

If you really want to dig into the chain and be sure it has fully converged, you can run the `auto_corrs` app over the MCMC chain:  

>> ./bin/auto_corrs output_tree.root

This will produce a root file, `output_tree_autoCorr.root` (if running over `output_tree.root`) containing a plot for each parameter of autocorrelation of parameter value at different lag (from 0 to 1000 steps).  

This can be useful step-size tuning and debugging.  

<h2> Physics Results </h2>

To be added!
 