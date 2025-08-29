<h1>Antinu Likelihood Analysis How To Guide</h1>

This file is a guide on how to run the most commonly used apps in the antinu oxo likelihood analysis.  

The first thing you'll want to do is get the environment variable in `env.sh` pointing to your `OXO` install. It's recommended you make a copy of this env.sh to avoid merge conflicts down the line.

Most of the executables live within the `exec` directory. Within `./src`, there is code for interfacing with config files and some other useful classes. Templates of the config files themselves live in `./cfg`.

<h2>Source Code and Classes</h2>

<h3>Config Loaders (src/config)</h3>

For each type of config file, there is a class defined in `./src/config`, along with a loader for that class. Each of these follows a similar structure: the config class has an attribute for each field in the config file, and the config loader classes read values in from the config files and return an instance of the corresponding config class.

<h3>Systematic Related Functions (src/syst)</h3>

<h4>Systematic Factory</h4>

The systematic factory, when given a systematic name, type, and list of associated fit parameters, will return a systematic object to be used in the analysis. For some systematics, a vector of oscillation grids and a reactor-distance map can also be handed to the systematic factory. The allowed types are the `OXO` systematic classes, along with a possible function name and ploy name, each separated with a colon. The possible functions and ploys are defined in the `SystematicFactory`. These are the current allowed combinations

- `Scale`: A linear Scale systematic
- `Shift`: A Shift systematic
- `Conv:Gaussian:SquareRootScale`: A Convolution with a Gaussian function where the width of the Gaussian runs with the square root of the axis in question
- `ScaleFunction:BirksLaw`: A ScaleFunction systematic where the function used modifies Birk's Constant
- `Shape:OscProb`: A Shape systematic, where the function used modifies the oscillation probability. The probability is calculated directly
- `Shape:OscProbGrid`: A Shape systematic, where the function used modifies the oscillation probability. The probability is obtained from OscGrids

<h4>Functions</h4>

Functions that can be passed to systematics are defined in here. Currently they are Birk's Law and the Oscillation Probability calculations. However, to pass lambdas to these functions, they had to be defined inline so this file may soon be removed.

<h4>OscGrids</h4>

The oscillation probability calculation can be quite slow to do on the fly, so instead we can produce an "Oscillation Grid" (<code>OscGrid</code>). The probability of oscillation depends on the true neutrino energy, the distance travelled, and the oscillation parameters $\Delta m^2_{21}$ and $\theta_{12}$, so we can pre-calculate the probabilities for different values of each of those variables, and the linearly interpolate between them. As the distances to reactor cores are all fixed distinct values, rather than interpolate, we can produce a 3D grid for each reactor core distance.

The `OscGrid` class has attributes for the min and max values and number of grid points for each dimension (true neutrino energy, $\Delta m^2_{21}$ and $\theta_{12}$), along with a vector of probabilities where each element represents the probability at a single point on the grid. An `OscGrid` can be constructed using just the axis ranges and number of points, along with the filename the grid will be written to and the distance it applies for. You can then calculate the probabilities at each grid point using the `CalcGrid` method. This can be written to the output file with the `Write` method. Alternatively you can load up a pre-written and calculated grid with the `Load` method. You can then obtain the probability for a given set of parameters with the `Evaluate` method. There are also 'Getters' for the probability vector, and vectors containing the grid points for each dimension. The calculation of the probability uses code copied from the RAT-tools/AntinuTools repo. We will look at making that code portable so we can call it rather than copy it in the future. The interpolation for calculating probabilities between grid points is done using functions which are stolen from the `OXO` Solar Analysis (https://github.com/dcookman/solar_analysis/). 

<h3>Utils (src/util)</h3>

<h4>Utilities</h4>

This file contains utility functions, mostly related to the Oscillation Grids. There is a function to turn the `OscGrids` max, min, and number of point values for a given dimension into a vector, along with a function to find the grid points that surround a given point in the `OscGrid`. These are both stolen from https://github.com/dcookman/solar_analysis. There is also a function here to load the reactor core indices and distances from SNO+ into a map from a JSON file, and functions taken from the rat-tools antinu tools to calculate the distance from reactors to SNO+.

<h4>DistBuilder</h4>

The DistBuilder class can build a <code>BinnedED</code> from a combination of a PDF Config and a dataset.

<h3>Configs (src/config)</h3>

There are several config files you'll need for running different apps. It's intended that you make a copy of the template config files and fill in the filepaths yourself. Just be aware of any subsequent changes to the format or content of the templates, as you'll probably want to propagate these to your own local config files.

<h4>Event</h4>

This file should contain information on the types of events you want to include in the analysis. An example is shown in `cfg/event_config.ini`. There should be a section, `summary`, which controls which event types are active and where the MC files live. Then there's a table for each event type, which contains a Latex name for axis titles, the exact location of files, the dimensions of the PDF, and the groups the PDF should belong to (different syestematics apply to different groups, see the `OXO` documentation for more info).

<h5>Summary</h5>

- `datasets`: The name of the datasets you want to fit. These should be comma separated, and each should be the name of a dataset table also defined in the file. 

<h5>Datasets</h5>

- `active`: The event types to be used in the analysis for this dataset. These should be comma separated, and each should also be the name of an Event field also defined in the file. You can also use `all` to use all events defined in the file
- `inactive`: Event types defined in the config that aren't being used for this dataset at this time
- `orig_base_dir`: The directory where the original (unpruned) MC files have been downloaded to
- `pruned_ntup_dir`: The directory where the pruned MC files should be written to, and later read from
- `datafile`: The file containing the data events for this dataset

<h5>Event Type Tables</h5>

The name of these tables should be the name of the event type.
- `ntup_files`: File path for the original unpruned files, relative to `orig_base_dir`
- `dimensions`: Number of dimensions the PDF should have
- `groups`: The groups the PDFs are part of. Each systematic is applied to one group (or all event types)

<h4>Fit</h4>

This file should contain information on the parameters of the fit, and the fit itself. An example is shown in `cfg/fit_config.ini`.

<h5>Summary</h5>

- `asimov`: Bool to determine if an Asimov dataset should be used. If 1/True, the scaled PDFs are used as the data to run an Asimov fit
- `fake_data`: Bool to determine if a 'fake dataset' should be used. If 1/True, PDFs are scaled to fake data values specified later in this config to produce the dataset. Fake data values of systematics are also applied. If both this and 'asimov' are true, the Asimov dataset is used
- `minuit_strategy` = The Minuit strategy level
- `minuit_tolerance` = The Minuit tolerance level
- `minuit_method` = The Minuit method used (`Migrad`, `Simplex` or `Minimize`)
- `iterations`: The number of steps in the Markov Chain/fit
- `burn_in`: The number of steps rejected at the start of the Markov Chain. This can be changed post analysis for most purposes
- `fit_dists`: Which PDFs should be fit. Can use `all`
- `output_directory`: Where to write output files to
- `n_steps`: Number of steps per iteration for Hamiltonian MCMC
- `epsilon`: Hamiltonian step size
- `hmc_iterations`: Number of steps in Hamiltonian MCMC
- `hmc_burn_in`: The number of steps rejected at the start of the Hamiltonian Markov Chain. This can be changed post analysis for most purposes
- `sigma_scale`: Global scaling factor each step size gets multiplied by
- `beeston_barlow`: Bool to determine if we should use the Beeston-Barlow method to account for MC stats uncertainty. If 1/True, it is used
- `save_outputs`: Bool to determine if the output files should be written. If 1/True, they are saved. This is useful for running multiple fixed oscillation fits to avoid producing ~millions of files at once. It is currently only implemented for the `fixedosc_fit` exec (for all other execs output files will be produced regardless of this bool)

<h5>Fit Parameter Tables</h5>

The name of these tables should be the name of the fit parameter.
- `nom`: Nominal value of the parameter
- `min`: The minimum value the parameter is allowed to take
- `max`: The maximum value the parameter is allowed to take
- `sig`: The relative width of the Gaussian used to propose new step values (each gets scaled by summary:sigma_scale). Often referred to as the step size
- `nbins`: Number of bins used for plotting the parameter
- `tex_label`: A latex name used for the parameter
- `constraint_mean`: Mean of the prior constraint on the parameter
- `constraint_sig`: Width of the prior constraint on the parameter
- `constraint_ratiomean`: Mean of prior constraint on ratio of this parameter to another
- `constraint_ratiosigma`: Width of prior constraint on ratio of this parameter to another
- `constraint_ratioparname`: Name of other parameter the ratio constraint corresponds to. This must be another parameter defined and activated in this file.
- `constraint_corr`: Correlation of this parameter to another
- `constraint_corrparname`: Name of other parameter this parameter is correlated to. This must be another parameter defined and activated in this file. Both correlated parameters must have a `constraint_mean` and `constraint_sig` defined.
- `fake_data_val`: Value to use to produce a fake dataset if `fake_data` is true in the Summary

A note on HMCMC: There is currently the functionality to run some MCMC, followed by some HMCMC starting from the best fit points of the MCMC. However, this didn't give any improvement in results or efficiency over the straight MCMC. The functionality is preserved in case we ever want to use it in the future, but generally we'll just run MCMC, with a notional 1 step HMCMC ran to avoid problems with files not being created when they are expected to.
It's likely in the future we'll just remove the functionality.

A note on Ratio Constraints: There is currently functionality to have a constraint on the ratio between two parameters. This should be specified by defining the constraint for one parameter and setting `constraint_ratioparname` to the name of the other. Do not then define the constraint for the second parameter as well, else the constraint will be applied twice. It does not matter which of the two parameters you define it for. It is not currently possible to have a ratio constraint for a parameter with more than one other parameters. We could change this in the future if it is needed.

Similarly, for correlations, do not define the correlation for the second parameter as well or it will be applied twice. It does not matter which of the two parameters you define the correlation for. It is not currently possible for a parameter to have a correlation with more than one other parameter. This could also change in the future.

<h4>PDF</h4>

This file should contain information on the pdf axes. An example is shown in `cfg/pdf_config.ini`.

<h5>Summary</h5>
  
- `build_order`: The axes used for the PDFs. These will be defined in the same file as axis tables. The order is important: the dimensions field in `events` config sets the number of dimensions, N, for a given PDF. The axes used will be the first N axes here
- `data_axes`: The axes the data will have

<h5>PDF Axis Tables</h5>
  
The name of these tables should be the name of the axis.
- `nbins`: The number of bins in the axis
- `min`: The minimum value on the axis
- `max`: The maximum value on the axis
- `branch_name`: The name of the branch in the ntuple
- `tex_name`: A latex name for the axis to be used in plot labels

<h4>Syst</h4>

This file should contain information on the systematics used in the analysis. An example is shown in `cfg/syst_config.ini`.

<h4>Summary</h4>
 
- `active`: Which systematics (defined in the systematics tables) will be used in the analysis

<h5>Systematic Tables</h5>

The name of these tables should be the name of the systematic.
- `param_names`: Comma separated list of the names of the fit parameters that control the systematic. These should all match a fit parameter defined in the fit config. If this is empty, the systematic name is assumed to be the name of the fit parameter
- `dist_obs`: These are the observables that the PDFs the systematic will apply to have
- `trans_obs`: These are the observables that the systematic will use to make the transformation
- `type`: The type of systematic being applied. This is not just the OXO derived systematic class, as here we also determine which function we will use for function-based systematics. The options are defined in the SystFactory (this will be made clearer in the future)
- `function`: The name of the function used if the systematic needs one. There is some redundancy here with type at the moment but we'll improve that soon
- `group`: This is the group of pdfs the systematic will be applied to. If this is left blank it will apply to all PDFs
- `dataset`: This is the dataset the systematic applies to. This should match the name of a dataset defined in the event config file

<h4>Oscillation Grid</h4>

This file should contain information about the oscillation grids you want to use in the analysis, both for producing them, and reading them.

<h5>Summary</h5>

- `filename`: Where the oscillation grids should be written to and read from
- `reactorsjson`: Where the reactor distances should be read from
- `mine`: The minimum true neutrino energy of the grid
- `maxe`: The maximum true neutrino energy of the grid
- `numvalse`: The number of grid points on the energy axis
- `mindm21sq`: The minimum $\Delta m^{2}_{21}$ value of the grid
- `maxdm21sq`: The maximum $\Delta m^{2}_{21}$ value of the grid
- `numvalsdm21sq`: The number of grid points on the $\Delta m^{2}_{21}$  axis
- `minssqth12`: The minimum $\text{sin}^2 \theta_{12}$ value of the grid
- `maxssqth12`: The maximum $\text{sin}^2 \theta_{12}$ value of the grid
- `numvalsssqth12`: The number of grid points on the $\text{sin}^2 \theta_{12}$ axis

<h2>Apps</h2>

These above classes and config files, along with all the OXO classes, are brought together in various apps in the `exec` directory. These are described in this section in the order you will probably want to use them:

<h3>Calculating Reactor Distances</h3>

To calculate the oscillation probability for an MC reactor IBD event, we need to know how far it has travelled, so we need to know the distance from the reactor it was produced in to SNO+. Rather than we do this calculation on the fly at every event for every iteration of a fit, we precalculate the distance for each reactor core. You can do this by running:

> ./bin/make_reactor_json osc_grid_config_file

This will loop over cores defined in the REACTORS RATDB table, and calculate the distance from it to SNO+. This gets saved in a JSON file, along with the reactor's name, and a unique integer assigned to the core. The integer is used as a numerical dimension for the reactor neutrino PDFs.

The output json file path is set in the oscillation grid config. These are likely not to change very often (only when cores come online, or an error in the locations in the RATDB table is fixed. If a reactor's distance to SNO+ physically changes we probably have more pressing issues!), so they can be saved and committed in the `reacttorjsons` directory to avoid having to reproduce them. It will be good practise to indicate in the filename the RAT version it was produced with, and which cores it contains. For RAT 8.0.0, there is a committed json file with all cores, and another one with just Bruce 1 for quicker testing of the machinery.

<h3>Making Oscillation Grids</h3>

To produce a single reactor core `OscGrid`, run:

> ./bin/make_osc_grids oscgrid_config index_to_write

In the submitting batch jobs section we will discuss running over all reactor cores at once.

<h3>Pruning Trees</h3>

Raw SNO+ ntuples, although much more lightweight than RATDS files, still contain much more information than we need at analysis level. So instead of carrying around that deadweight for the whole analysis, we first take the time to prune out the branches we don't need and save what we do need in new files to be used.  

The `prune_trees` app does this for us. Having set your filepaths for each dataset in your `events` config file, simply run:  

> ./bin/prune_trees cfg/event_config.ini

This looks for ntuple files for each of the event types specified in the `event` config (you can select/deselect event types with the active/inactive fields) and produces a new set of files containing only the relevant trees for the analysis (energy, position, fit validity, pulse shape discrimination values, and true neutrino energy and reactor name for reactor neutrinos). Be warned, if using all event types, this will take some time! (of order several hours).

This app also divides the three alpha-n processes into different files. They are all simulated at once, and so each of the three alpha-n processes gets events from the same files. For each alpha-n process, we loop over all alpha-n events, and select the events of that specific alpha-n process using the reconstructed energy.

<h3>LLH Scans</h3>

Right, we've now pruned our trees and can use them to build PDFs and a simulated dataset(s). In theory, we're now ready to launch a fit! But let's be steady Eddies, Cautious Carols, and Nervous Nigels, and make sure everything is behaving as intended. Running fits can be computationally and storage expensive, so it would be a shame to find out after running them that something had previously gone wrong.  

One way of doing this is with likelihood scans. The likelihood scan apps first build the same test statistic that the fit will use, and we will hand it the same PDFs, parameters, and systematics we intend to hand the fitter. For a single floated parameter, it will scan over a range of values centred at the nominal value. For each of 150 steps, it will rebuild the dataset by scaling the PDFs and applying any systematics all set at their nominal values, except the parameter in question which is set to the value we've reached in the scan. This is compared to the actual dataset using the test statistic (probably the binned LLH). The LLH (minus the nominal LLH) is saved at each step, and once we reach the end of the scan for a parameter, it is set back to its nominal value, and we repeat for the next parameter.  

For the 'true' Asimov dataset, these scans should all minimise at (1,0) (where the x axis is parameter value / nominal value) i.e., by changing a parameter, you can't produce a dataset more similar to the target dataset than by using the exact values used to produce the target dataset. For the other forms of Asimov dataset, most parameters will minimise close to 1, but probably not exactly. The event types with higher rates will be closer to 1 as changes in these will have a greater impact on the LLH, and the small fluctuations causing the differences between the PDFs and Asimov dataset will be smaller proportionally.  

There are two likelihood scan apps. The first, `fixedosc_llhscan` initially loops over all parameters except the oscillation parameters. The oscillation parameters are not handed to the `OXO` `BinnedNLLH` object, and the reactor PDF is produced with the oscillation parameters at nominal values. After performing the LLH scan for all non-oscillation parameters, it rebuilds the reactor PDF for each of the 150 values for each oscillation parameter. For each of these, it rebuilds the dataset and hands this to the `BinnedNLLH` to compare to the Asimov dataset. This is because, for speed, sometimes we want to run fits where the oscillation is handled outside of `OXO` (see the fixed oscillation parameters fits below). This version of the likelihood scan is designed to be used for validating that fit.

To run the fixedosc LLH scan:  

> ./bin/fixedosc_llhscan cfg/fit_config.ini cfg/event_config.ini cfg/pdf_config.ini cfg/syst_config.ini cfg/oscgrid.ini

You will need to have a `deltam21` parameter, and exactly one `theta12`, `sintheta12`, or `sinsqtheta12` parameter specified in the fit config. A root file `llh_scan.root` will be saved in the output directory specified in the fit config. In the file will be a plot of LLH vs parameter value for each parameter.

The other version of the likelihood scan, `llh_scan`, hands all the parameters, including the oscillation parameters, to the `BinnedNLLH`, but it is currently deprecated (and not compiled by default). It scans through all parameters in the same way, setting the values in the `BinnedNLLH` and evaluating the likelihood. This version of the likelihood scan is designed to be used for validating the full fits where every parameter floats at once (see below).

To run the LLH scan:  

> ./bin/llh_scan cfg/fit_config.ini cfg/event_config.ini cfg/pdf_config.ini cfg/syst_config.ini cfg/oscgrid.ini

Both LLH scan apps will also produce several subdirectories inside the `output_directory` set in the `fit_config`. `unscaled_pdfs` will contain each of the unscaled (normalised) PDFs without any systematics applied. `asimov_dists` will contain each of the PDFs scaled to their nominal rates with nominal systematics applied. If you're running with `fake_data = 1` in the `fit_config`, `fakedata_dists` will contain each of the PDFs scaled to their fake data value rate with systematics applied at their fake data values.

<h3>Running a Fit</h3>

Now the time has come to run a fit. When we float the oscillation parameters, we need extra dimensions for the reactors PDF to contain the true neutrino energy and the reactor core it came from. As there are 483 cores, if we have 140 bins in each energy dimension, we end up with over 9 million bins! This means marginalising down to 1D to compare to data takes a lot of time. There are few ideas to get around this, but first let's run individual fits where we don't vary the oscillation parameters (so everything remains 1D and fast), and do this at all points in a grid of oscillation parameter values. We can do this for one point with:

> ./bin/fixedosc_fit cfg/fit_config.ini cfg/event_config.ini cfg/pdf_config.ini cfg/syst_config.ini cfg/oscgrid.ini

The value of the oscillation parameters should be set in the `fit` config file. You will need to have a `deltam21` parameter, and exactly one `theta12`, `sintheta12`, or `sinsqtheta12` parameter. This performs a `Minuit` fit, with the oscillation parameters fixed at those values and the Minuit settings in the `fitCofig`. It will load up all the pdfs, systematics, data etc. and creates the likelihood object in the same way the `llh_scan` app does, and then runs a fit. A variety of output files are produced (if the `save_outputs` in the `fitConfig` is set).

`fit_results.txt` contains each parameter value for the maximum LLH point, and `fit_results.root` contains the vectors of parameter names, postfit values and uncertainties, and a covariance matrix.

There will also be several subdirectories produced inside the `output_directory` set in the `fit_config`. `unscaled_pdfs` will contain each of the unscaled (normalised) PDFs without any systematics applied. `asimov_dists` will contain each of the PDFs scaled to their nominal rates with nominal systematics applied. If you're running with `fake_data = 1` in the `fit_config`, `fakedata_dists` will contain each of the PDFs scaled to their fake data value rate with systematics applied at their fake data values.

`postfit_dists` contains the PDFs for each event type, scaled by the event type parameter value at the maximum LLH step in this chain, with each systematic’s value at the maximum LLH step applied. Also saved is the sum of these.

In `1dlhproj` and `2dlhproj`, you have projections of the LLH for each parameter, and each combination of two parameters. In each case, all other parameters are marginalised over.

If you want to run a full fit where everything is floated at once, there is an app to run some Markov Chain Monte Carlo. First you want to make sure everything in your fit config is looking sensible. The min and max values, and sigmas for each event type are fairly well tuned, so it's advised you leave these alone unless you know what you're doing. Certainly, for your first fit testing out the code I would leave it as is. If you have good reason to change them, you probably don't need to be reading this guide! You can run a full fit by:

> ./bin/mcmc cfg/fit_config.ini cfg/event_config.ini cfg/pdf_config.ini cfg/syst_config.ini cfg/oscgrid.ini

This will load up all the pdfs, systematics, data etc. and creates the likelihood object in the same way the `llh_scan` app does, and then runs the MCMC. For each individual chain, you’ll have a number of files and subdirectories, some of which are the same as for `fixedosc_fit`.

`auto_correlations.txt` contains how correlated the LLH is to the LLH at a step a different number of steps previous. This should hopefully be close to 0 after a few thousand steps (there is a separate `autocorrelations` app you can use to run over the outputted tree and get more information). `fit_results.txt` contains each parameter value for the maximum LLH step.

There will also be several subdirectories produced inside the `output_directory` set in the `fit_config`. `unscaled_pdfs` will contain each of the unscaled (normalised) PDFs without any systematics applied. `asimov_dists` will contain each of the PDFs scaled to their nominal rates with nominal systematics applied. If you're running with `fake_data = 1` in the `fit_config`, `fakedata_dists` will contain each of the PDFs scaled to their fake data value rate with systematics applied at their fake data values.

`postfit_dists` contains the PDFs for each event type, scaled by the event type parameter value at the maximum LLH step in this chain, with each systematic’s value at the maximum LLH step applied. Also saved is the sum of these.

In `1dlhproj` and `2dlhproj`, you have projections of the LLH for each parameter, and each combination of two parameters. In each case all other parameters are marginalised over. The burn-in steps are automatically not included.

There will also be a root file (`fit_name_i.root`) which contains a tree. Each entry in the tree represents a step in the Markov Chain, and the leaves are the values of each parameter, as well as the LLH, step time, and acceptance rate. This is more for MCMC chain diagnostics and debugging than obtaining physics results but is useful for checking the fit has worked as expected.

<h2>Submitting Batch Jobs</h2>

<h3>Individual Jobs</h3>

Every batch system will be different, but for submitting to a Condor based queue system there is a script, `util/submitCondor.py`, based on one for submitting RAT jobs originally from Josie (I think). You can use this for any exectuable. For running full MCMC fits, it’s advisable to run multiple identical fits at once in parallel as you probably want around 1 million steps.   

You can submit N jobs with:  

> python utils/submitCondor.py exec_name output_dir -r /path/to/this/repo/ -e /path/to/env/file/ -f cfg/fit_config.ini -i cfg/event_config.ini  -p cfg/pdf_config.ini -s cfg/syst_config.ini -o cfg/oscgrid.ini -n N -w walltime

The environment file you supply here should be whatever you use to set up ROOT, python, GSL, etc., not `antinullh/env.sh` (which should be sourced before running the python submission). The output directory will be created, with several subdirectories:

`error`, `output` and `log` contain the console outputs and batch log. `sh` and `submit` contain the files used to submit these jobs, and `cfg` contains all the config files used to run the fit. The `output_directory` set in the `fit_config` file will be overwritten by the command line `output_dir` argument.

You can run this script with any of the apps. If you're not running a fit, you probably don't need multiple simultaneous jobs so can just run with N=1. If you're running `make_osc_grid` with the submission script, it will automatically loop over every reactor core in the reactors JSON file and run `make_osc_grid` for each in parallel.

<h3>Fixed Oscillation Parameters Fits</h3>

First, you may want to do a grid scan of oscillation parameters, running a fit at each step. In these fits, the oscillation parameters are fixed at the values for that point in the scan, and everything else is floated in a Minuit fit. We would normally do ~500 steps for each oscillation parameter, so 250,000 individual fits. Now, you can try submitting 250,000 jobs at once but do so at your own (and your friendly neighbourhood sys-admin's) peril! Instead, we do one job for each value of one of the oscillation parameters, so that's 500 jobs each running 500 sequential fits. This number can be tuned for efficient running on your own cluster. The -d option allows you to change the number of fits ran per job.

There is a script, similar to `utils/submitCondor.py` but for submitting these fixed oscillation parameter fits. This is in `util/submitFixedOscJobs.py`. 

You can submit the jobs with:

> python utils/submitFixedOscJobs.py fixedosc_fit output_dir -r /path/to/this/repo/ -e /path/to/env/file/ -f cfg/fit_config.ini -i cfg/event_config.ini  -p cfg/pdf_config.ini -s cfg/syst_config.ini -o cfg/oscgrid.ini -d fits_per_job -w walltime (s) -m requested_memory (MB) -j jobname

In your fit config, be sure to have `fake_data=1` (and `asimov=0`), and set the oscillation parameter's fake data values to be their nominal values. This way, for every step in the grid, the fake dataset will be produced with the nominal values, and this will be the same for each fit. The python script will update the nominal values of the oscillation parameters in the config file to the current values in the fixed oscillation parameters scan, so we always fit the same fake data, but with different fixed oscillation parameter values each time. 

`output_dir` will contain all the configs and logs for each fit, and the outputted root files will be produced in independent directories (one for each fit) inside that directory. The `output_directory` set in the `fit_config` file gets updated to be the one set as the command line argument.
 
<h2>Postfit Analysis</h2>

<h3>Fixed Oscillation Fits</h3>

A lot of the postfit analysis will be done using scripts that reside in `util` (not to be confused with `src/util`), and `plotting`.

<h4>makeFixedOscTree</h4>

The first thing you'll want to do after running a set of fixed oscillation parameter fits is run `exec/makeFixedOscTree.cc`. This is currently the only postfit analysis step that requires compiling (which will be handled by the usual Makefile). It loops over the outputs of all the fits, and fills a tree. Each entry in the tree represents one fit, and each branch is a fit parameter. The nominal values and prefit constraints are also saved in vectors. You can run it with:

> ./bin/makeFixedOscTree cfg/fit_config.ini cfg/oscgrid_config.ini

The config files should be ones you've used to run one of the fits. If this takes a long time, you can submit it as a batch job using `utils/submitCondor.py`:

> python utils/submitCondor.py makeFixedOscTree output_dir -r /path/to/this/repo/ -e /path/to/env/file/ -f cfg/fit_config.ini-o cfg/oscgrid.ini -w walltime

Once the tree is made, the it finds the fit that had the best LLH, and reruns that fit with the `save_outputs` bool in the fit config set to true and Minuit settings updated to do a more rigorous fit. The covariance matrix from this fit is saved along with the tree of fit results.

<h4>plotFixedOscLLH</h4>

The next thing to do is to plot the best fit LLH as a function of the oscillation parameters. You can do this by running `plotting/plotFixedOscLLH.C` over the output `TTree` from `makeFixedOscTree`. It loops over all the entries, and gets the oscillation parameter values and LLH for each fit. It then plots a 2D histogram where the X and Y axis are the oscilltion parameters, and the Z axis is the LLH. A canvas is saved in both a `.root` and `.pdf` file, in the top level output directory of the set of fits. 

Also plotted are profile LLHs for both $\Delta m^2_{21}$ and $\theta_{12}$. Each are saved as `.pdf` file and also in the same root file as the 2D plot. The $1\sigma$ bounds are calculated and saved in that root file. These can then be used by `plotFixedOscParams` as the oscillation parameter postfit uncertainties.

It can be run by doing:

> root -l 'plotting/plotFixedOscLLH.C("/path/to/makeFixedOscTree/output")'

<h4>plotFixedOscDist</h4>

This script will loop over all entries in the output `TTree` from `makeFixedOscTree`, and find the fit with the minimum best LLH. It then goes to the directory of that fit, and plots the distributions saved in the `postfit_dists` directory. The distributions plotted are the data, the total MC (sum of all scaled PDFs with systematics applied), and each individual PDF scaled (without systematics applied). Also saved is a similar plot but with PDFs grouped together. For both plots, a panel showing the ratio of the total MC postfit prediction to the data is plotted below the main histogram.

Currently, the user hard codes the parameter names for each dataset in the top of this file. When running, the user supplies an option as the second argument. `1` plots the distributions for the first dataset, `2` plots the distributions for the second dataset, and `0` plots the summed distributions for both datasets. Obviously it would be nice if this wasn't hardcoded, and you could use N datasets, but for expedience this is how it is at the moment.

Canvases are saved in both `.root` and `.pdf` files, in the top level output directory of the set of fits. You can run it with:

> root -l 'plotting/plotFixedOscDist.C("/path/to/makeFixedOscTree/output", dataset_option)'

In this script Latex labels are made for each PDF (currently reactor IBDs, Geo U, Geo Th, the three #alpha-ns, and Sideband). If more PDFs are used, these will be added to the backgrounds group without a Latex label. It is recommended this script gets updated if the fit parameters change.

<h4>compare3FixedOscDists</h4>

This script is similar to `plotFixedOscDist`, but plots the total postfit scaled MC for 3 fits (without each scaled component PDF). It loops through the entries in the output `TTrees` from `makeFixedOscTree`, and finds the fit with the minimum best LLH for each. The distributions plotted are the data for the first fit, and the total MC (sum of all scaled PDFs with systematics applied) for each of the three fits, and a panel showing the ratio of each total MC postfit prediction to the data for the first fit.

> root -l 'plotting/compare3FixedOscDists.C("/path/to/fit1/params.root", "/path/to/fit2/params.root", "/path/to/fit3/params.root", "fit 1 label", "fit 2 label", "fit 3 label")'

<h4>plotFixedOscParams</h4>

This script loops over all entries in the output `TTree` from `makeFixedOscTree`, and finds the fit with the minimum best LLH. It then plots each parameter in that entry, relative to it's nominal value. Any prefit constraints are also plotted, relative to nominal values. A canvas is saved in both a `.root` and `.pdf` file, in the top level output directory of the set of fits. The constraint, nominal, and postfit histograms are also saved in the root file. If `plotFixedOscLLH` has already been run, it will look for the outputted `LLH.root` file to find the oscillation parameter uncertainties. Otherwise, it will use the grid spacing.

The reactor parameter names are hardcoded so that they can be multiplied by the correct reactor ratio. The order the parameters are plotted in are also hardcoded.

A separate canvas is also saved, containing a plot of the postfit correlation matrix of all fit parameters.

You can run it with:

> root -l 'plotting/plotFixedOscParams.C("/path/to/makeFixedOscTree/output")'

<h4>compareFixedOscParams</h4>

This script opens the outputted root file from `plotFixedOscParams` for two different fits, and plots the postfit parameters for each and the prefit nominals and constraints from the first. The user also inputs labels for each fit to be printed into the legend, along with the name (without file type suffix) of the outputted files. You can run it with:

> root -l 'plotting/compareFixedOscParams.C("/path/to/fit1/params.root", "fit 1 label", "/path/to/fit2/params.root", "fit 2 label", "/path/to/output/file")'

<h4>compare3FixedOscParams</h4>

This is just like `compareFixedOscParams` but it runs over 3 fits:

> root -l 'plotting/compare3FixedOscParams.C("/path/to/fit1/params.root", "fit 1 label", "/path/to/fit2/params.root", "fit 2 label", "/path/to/fit3/params.root", "fit 3 label", "/path/to/output/file")'. 

<h4>compare2LLHScans</h4>

This script loops through the objects in the first of two LLH scan output files and prints each histogram to a canvas. If a histogram with the same name exists in the second LLH scan output files, these are drawn on the same canvas. Each canvas is saved as one page in a PDF file, and in a root file. You can run it with:

> root -l 'plotting/compare2LLHScans.C("/path/to/llhscan1/llhscan.root", "/path/to/llhscan1/llhscan.root", "fit 1 label", "fit 2 label")'

<h4>compare3LLHScans</h4>

This script is the same as `compare2LLHScans` but for three scans. It loops through the objects in the first LLH scan output file and prints each histogram to a canvas. If a histogram with the same name exists in the second and/or third LLH scan output files, these are drawn on the same canvas. Each canvas is saved as one page in a PDF file, and in a root file. You can run it with:

> root -l 'plotting/compare3LLHScans.C("/path/to/llhscan1/llhscan.root", "/path/to/llhscan1/llhscan.root", "/path/to/llhscan3/llhscan.root", "fit 1 label", "fit 2 label", "fit 3 label")'

<h4>plot1DPDFs</h4>

This script loops through all the unscaled PDF files in a directory, and plots each on a different page of a PDF file, and in a root file. You can run it with:

> root -l 'plotting/plot1DPDFs("/path/to/dir/")'

<h4>compare1DPDFs</h4>

This script loops through all the unscaled PDF files in a directory, and plots each on a different page of a PDF file, and in a root file. The user inputs the two directory paths, and labels for the legend. You can run it with:

> root -l 'plotting/compare1DPDFs("/path/to/dir1/", "/path/to/dir2/", "label1", "label2")'

<h4>plot2DPDFs</h4>

This script loops through all the 2D unscaled PDF files in a directory, and plots each (with colz) on a different page of a PDF file, and in a root file. You can run it with:

> root -l 'plotting/plot2DPDFs("/path/to/dir/")'

<h4>projectPDFs</h4>

This script loops through all the 2D unscaled PDF files in a directory, and plots the projection of each on the y-axis (likely to be the #alpha-n classifier axis) on a different page of a PDF file, and in a root file. You can run it with:

> root -l 'plotting/projectPDFs("/path/to/dir/")'

<h4>compareContours</h4>

This script uses the files created by `plotFixedOscLLH` and draws the 1 $\sigma$ contours for two fits. The user inputs the two file paths and labels for the legend, along with the name (without file type suffix) of the outputted files.

> root -l 'plotting/compareContours("/path/to/fitdir1/plots/LLH.root", "/path/to/fitdir2/plots/LLH.root", "label1", "label2", "/path/to/output/file")'

<h3>MCMC</h3>

NOTE: This section is a bit out of date because the code hasn't been brought in line with the rest of the antinullh repo yet. It will be updated when the code is! The non-fixed osc fits and llh scans are deprecated for now, and are not compiled by default.

Now you’re fit has run, the real fun starts! For each individual chain, you’ll have a number of files and subdirectories. In `1dlhproj` and `2dlhproj`, you have projections of the LLH for each parameter, and each combination of two parameters. In each case all other parameters are marginalised over. The burn-in steps, as set in the fit config, are automatically not included.  

`config_log.txt` and `config_log_hmc.txt` save the configs used (these are almost certainly the same as each other). These are also saved in the `cfg` directory if running the python submission script in `util` (see below). `auto_correlations.txt` calculates how correlated the LLH is to the LLH at a step a different number of steps previous. This should hopefully be close to 0 after a few hundred. `fit_results.txt` contains each parameter value for the maximum LLH step. `postfit_dists` contains the PDFs for each event type, scaled by the event type parameter value at the maximum LLH step in this chain. Also saved is the sum of these, with each systematic’s value at the maximum LLH step applied.  

There will also be a root file (`fit_name_i.root`) which contains a tree. Each entry in the tree represents a step in the Markov Chain, and the leaves are the values of each parameter, as well as the LLH, step time, and acceptance rate. This is more for MCMC chain diagnostics and debugging than obtaining physics results but is useful for checking the fit has worked as expected.  

To combine multiple parallel chains, `hadd` these root files. When plotting from the `hadded` file, you can cut off the burn-in using the `Step` branch (as entry will now run from 0 to total number of combined steps, not 0 to number of steps in each file).  

There’s an app for making useful plots from the combined tree. This can be run as:  

> ./bin/make_plots output_tree.root asimov_dataset.root scaled_postfit_dist.root fit_config_file.cfg correlations traces

`output_tree.root` is either a single outputted file or the `hadded` combination. `asimov_dataset.root` is the ‘true’ Asimov dataset to get the Asimov rates to compare to. `scaled_postfit_dist.root` is the sum of the PDFs scaled by the maximum LLH step values. If you’ve `hadded` chains, you can find the chain with the highest LLH step (so you know which `postfit_dists` to use) using `util/findMaxStep.C`:  

> root -l -b -q 'findMaxStep.C("dataset_directory_name", "fit_name", number_of_chains )'

`correlations` and `traces` determine if you draw the correlations between all combinations of two parameters (parameter 1 vs parameter 2 for all steps) and the trace of each parameter (parameter vs Step). These are both very time consuming. Leaving them empty means neither will run (note that the code just looks for a non-NULL string, so setting them to ‘0’ means they will be run!).  

The `make_plots` code is pretty messy, with features having been tacked on as time has gone on. But it does a job. The outputs are a root file of histograms, and a PDF with a bunch of histograms drawn to separate pages. For each parameter, you’ll get a LLH distribution marginalised over all other parameters. This is essentially a frequency histogram of parameter value for all steps. The central postfit value is calculated from this histogram in four different ways: the arithmetic mean of the histogram, the mean of a Gaussian fitted to it, the highest posterior density (mode), and the value at the maximum LLH step. These central values (and associated uncertainties) are combined into a single plot after the individual parameter plots. Prefit and postfit (using the maximum LLH step values) distributions are plotted and compared in one and two dimensions. If `traces` was set, a plot of parameter value vs `Step` is also drawn. If `correlations` are a run, a `colz` plot of parameter values at all steps are drawn.  

If you really want to dig into the chain and be sure it has fully converged, you can run the `auto_corrs` app over the MCMC chain:  

> ./bin/auto_corrs output_tree.root

This will produce a root file, `output_tree_autoCorr.root` (if running over `output_tree.root`) containing a plot for each parameter of autocorrelation of parameter value at different lag (from 0 to 1000 steps).  

This can be useful step-size tuning and debugging.  

<h2>Physics Results</h2>

To be added!
 
