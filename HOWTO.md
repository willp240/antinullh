<h1>Antinu Likelihood Analysis How To Guide</h1>

This file is a guide on how to run the most commonly used apps in the antinu oxo likelihood analysis.  

The first thing you'll want to do is set `OXO_DIR` and the other build environment variables. You can copy `env.sh` to a local file and adjust it without changing the shared version.

Most of the executables live within the `exec` directory. Within `./src`, there is code for interfacing with config files and some other useful classes. Templates of the config files themselves live in `./cfg`.

<h2>Source Code and Classes</h2>

<h3>Config Loaders (src/config)</h3>

For each type of config file, there is a class defined in `./src/config`, along with a loader for that class. Each of these follows a similar structure: the config class has an attribute for each field in the config file, and the config loader classes read values in from the config files and return an instance of the corresponding config class.

<h3>Systematic Related Functions (src/syst)</h3>

<h4>Systematic Factory</h4>

`SystFactory` creates one OXO systematic from its name, type, parameter names and nominal values. `BuildSysts` reads the syst config and creates the systematics needed for each dataset. Types can include a function and ploy name separated by colons. The supported types include

- `Scale`: A linear Scale systematic
- `Shift`: A Shift systematic
- `Conv:Gaussian:SquareRootScale`: A Convolution with a Gaussian function where the width of the Gaussian runs with the square root of the axis in question
- `ScaleFunction:BirksLaw`: A ScaleFunction systematic where the function used modifies Birk's Constant
- `Shape:ReactorOscProb`: A Shape systematic that applies reactor oscillation weights during the fit

`SystFactory.cc` also defines the available classifier and normalisation Shape functions. See its `Shape` cases for their type names.

<h4>Functions</h4>

`SystFuncs.hh/.cc` contains the reactor oscillation calculation and the function that prepares its Shape from binned MC. Some other Shape functions are defined where `SystFactory` constructs them.

<h3>Utilities</h3>

Reactor JSON loading and baseline calculations are in `src/config/ReactorBaselines.hh/.cc`. The two configuration string helpers are in `src/util/StringUtils.hh/.cc`.

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
- `data_file`: The file containing the data events for this dataset

<h5>Event Type Tables</h5>

The name of these tables should be the name of the event type.
- `ntup_files`: File path for the original unpruned files, relative to `orig_base_dir`
- `dimensions`: Number of dimensions the PDF should have
- `groups`: The groups the PDFs are part of. Each systematic is applied to one group (or all event types)
- `oscillated`: Should the events be oscillated before making the PDF. Defaults to false
- `flat`: Should the PDF ignore the MC events and make a flat uniform distribution

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
- `livetime`: Exposure setting read by the config loader; the current fit executables do not use it when building distributions
- `n_steps`, `epsilon`, `hmc_iterations`, `hmc_burn_in`: Legacy Hamiltonian MCMC settings still read by the config loader but not used by the current `mcmc` executable
- `sigma_scale`: Global scaling factor each step size gets multiplied by
- `beeston_barlow`: Bool to determine if we should use the Beeston-Barlow method to account for MC stats uncertainty. If 1/True, it is used
- `save_outputs`: Whether to save input spectra and detailed fit outputs. Scan executables always save their `llh_scan.root` result.
- `llh_scan_points`: For the llh scan execs, this will set how many points are in the scan (defaults to 150). Other execs will just ignore this parameter

<h5>Fit Parameter Tables</h5>

The name of these tables should be the name of the fit parameter.
- `nom`: Nominal value of the parameter
- `min`: The minimum value the parameter is allowed to take
- `max`: The maximum value the parameter is allowed to take
- `sig`: The relative width of the Gaussian used to propose new step values (each gets scaled by summary:sigma_scale). Often referred to as the step size
- `nbins`: Number of bins used for plotting the parameter
- `tex_label`: A latex name used for the parameter
- `constraint_mean`: Mean of the prior constraint on the parameter
- `constraint_sigma`: Width of the prior constraint on the parameter
- `constraint_ratiomean`: Mean of prior constraint on ratio of this parameter to another
- `constraint_ratiosigma`: Width of prior constraint on ratio of this parameter to another
- `constraint_ratioparname`: Name of other parameter the ratio constraint corresponds to. This must be another parameter defined and activated in this file.
- `constraint_fracmean`, `constraint_fracsigma`, `constraint_fracparname`: Mean, width and other parameter for a fraction constraint
- `constraint_shapemean`, `constraint_shapesigma`, `constraint_shapeparnames`, `constraint_shapefunc`: Settings for a shape constraint
- `constraint_corr`: Correlation of this parameter to another
- `constraint_corrparname`: Name of other parameter this parameter is correlated to. This must be another parameter defined and activated in this file. Both correlated parameters must have a `constraint_mean` and `constraint_sigma` defined.
- `fake_data_val`: Value to use to produce a fake dataset if `fake_data` is true in the Summary
- `fix`: Whether or not parameter is fixed in the fit

The current `mcmc` executable uses a Metropolis sampler. It does not run an HMCMC stage.

A note on Ratio Constraints: There is currently functionality to have a constraint on the ratio between two parameters. This should be specified by defining the constraint for one parameter and setting `constraint_ratioparname` to the name of the other. Do not then define the constraint for the second parameter as well, else the constraint will be applied twice. It does not matter which of the two parameters you define it for. It is not currently possible to have a ratio constraint for a parameter with more than one other parameters. We could change this in the future if it is needed.

Similarly, for correlations, do not define the correlation for the second parameter as well or it will be applied twice. It does not matter which of the two parameters you define the correlation for. It is not currently possible for a parameter to have a correlation with more than one other parameter. This could also change in the future.

<h4>PDF</h4>

This file should contain information on the pdf axes. An example is shown in `cfg/pdf_config.ini`.

<h5>Summary</h5>
  
- `build_order`: The axes used for the PDFs. These will be defined in the same file as axis tables. The order is important: the dimensions field in `events` config sets the number of dimensions, N, for a given PDF. The axes used will be the first N axes here
- `data_axes`: The axes the data will have

<h5>PDF Axis Tables</h5>
  
The name of these tables should be the name of the axis.
- `n_bins`: The number of bins in the axis
- `min`: The minimum value on the axis
- `max`: The maximum value on the axis
- `branch_name`: The name of the branch in the ntuple
- `tex_name`: A latex name for the axis to be used in plot labels
- `llh_buffer_bin_low`: Optionally set the last bin of the lower buffer region (defaults to 8 for energy, won't be set for other axes so no buffer is applied)
- `llh_buffer_bin_high`: Optionally set the first bin of the upper buffer region (defaults to 20 for energy, won't be set for other axes so no buffer is applied)

<h4>Syst</h4>

This file should contain information on the systematics used in the analysis. An example is shown in `cfg/syst_config.ini`.

<h4>Summary</h4>
 
- `active`: Which systematics (defined in the systematics tables) will be used in the analysis

<h5>Systematic Tables</h5>

The name of these tables should be the name of the systematic.
- `param_names`: Comma separated list of the names of the fit parameters that control the systematic. These should all match a fit parameter defined in the fit config. If this is empty, the systematic name is assumed to be the name of the fit parameter
- `dist_obs`: These are the observables that the PDFs the systematic will apply to have
- `trans_obs`: These are the observables that the systematic will use to make the transformation
- `type`: The OXO systematic class, optionally followed by a function and ploy name separated by colons, such as `Conv:Gaussian:SquareRootScale` or `Shape:ReactorOscProb`
- `group`: This is the group of pdfs the systematic will be applied to. If this is left blank it will apply to all PDFs
- `dataset`: This is the dataset the systematic applies to. This should match the name of a dataset defined in the event config file

<h4>Reactor Baselines</h4>

The fit config's `[summary]` section includes:

- `reactorsjson`: path to the reactor-index, name and baseline JSON. Relative paths are resolved from the executable's working directory.

<h4>Fixed-Oscillation Scans</h4>

The fit config's `[summary]` section also includes:

- `fixedosc_dm_points`: number of fixed-fit scan points in delta-m-squared (default 500).
- `fixedosc_theta_points`: number of fixed-fit scan points in the configured theta-like parameter (default 500).

The scan bounds come from the oscillation parameters' `min` and `max` fields. These scan counts describe repeated fixed fits; MCMC does not use them. Reactor oscillation probabilities are calculated directly, with no probability lookup tables.

<h3>Fit Code Layout</h3>
All four fit and scan executables live in `exec/`. They read the same four INI files and follow the same stages: load and check configuration, prepare dataset spectra, build the likelihood, run the fit or scan, and write results.

<h4>`src/config/`</h4>
- AnalysisConfig.hh/.cc – Loads the four INI files together and checks references between datasets, components, systematics, and fit parameters.
- ReactorBaselines.hh/.cc – Read the reactor baseline map and calculate distances when generating it.
- The config and loader classes – Parse each file's own fields.

<h4>`src/data/`</h4>
- DatasetSetup.hh and LoadEvents.hh – Describe datasets and load data and MC files.
- BuildDist.hh/.cc – Fill binned distributions from events, including the MC used for reactor oscillation.
- BuildComponent.hh/.cc – Apply systematics and rate scaling to a component PDF, then project it onto the data observables.
- BuildDatasets.hh/.cc – Assemble component PDFs and observed, Asimov, or fake data spectra. It does not write files or construct likelihood objects.

<h4>`src/fit/`</h4>
- FitParameters.hh – Working fit parameters, parameter selection and rate scaling.
- OscParams.hh – Interpret the configured oscillation parameters.
- BuildLikelihood.hh/.cc and CombineLLHs.hh – Build per-dataset OXO likelihoods and combine them with shared parameters.
- FitConstraints.hh and ShapeConstraints.hh – Apply parameter constraints.

<h4>`src/syst/`</h4>
- BuildSysts.hh/.cc – Construct and own the systematics for each dataset.
- ApplySysts.hh – Apply relevant systematics to a component by group.
- SystFactory.hh/.cc and SystFuncs.hh/.cc – Implement the systematic types and functions.

<h4>`src/output/`</h4>
- WriteDists.hh – Write input and fitted distributions
- OutputDirs.hh, FitSummary.hh and RateReporting.hh – Output paths and summaries.

<h2>Apps</h2>

These above classes and config files, along with all the OXO classes, are brought together in various apps in the `exec` directory. These are described in this section in the order you will probably want to use them:

<h3>Calculating Reactor Distances</h3>

To calculate the oscillation probability for an MC reactor IBD event, we need to know how far it has travelled, so we need to know the distance from the reactor it was produced in to SNO+. Rather than do this calculation on the fly at every event for every iteration of a fit, we precalculate the distance for each reactor core. You can do this by running:

> ./bin/make_reactor_json reactorjsons/allcores_rat_VERSION.json

This will loop over cores defined in the `REACTORS` RATDB table, and calculate the distance from it to SNO+. This gets saved in a JSON file, along with the reactor's name, and a unique integer assigned to the core. The integer is used as a numerical dimension for the reactor neutrino PDFs.

The output JSON path is supplied on the command line. Set `reactorsjson` in the fit config to that path. Keep the same JSON mapping for pruning and fitting: the stored reactor indices must correspond to the MC's `reactorIndex` values. These are likely not to change very often (only when cores come online, or an error in the locations in the RATDB table is fixed. If a reactor's distance to SNO+ physically changes we probably have more pressing issues!), so they can be saved and committed in the `reactorjsons` directory to avoid having to reproduce them. It will be good practice to indicate in the filename the RAT version it was produced with, and which cores it contains. For RAT 8.0.0, there is a committed json file with all cores, and another one with just Bruce 1 for quicker testing of the machinery.

<h3>Pruning Trees</h3>

Raw SNO+ ntuples, although much more lightweight than RATDS files, still contain much more information than we need at analysis level. So instead of carrying around that deadweight for the whole analysis, we first take the time to prune out the branches we don't need and save what we do need in new files to be used.  

The `prune_trees` app does this for us. Having set your filepaths for each dataset in your `events` config file, simply run:  

> ./bin/prune_trees cfg/event_config.ini cfg/fit_config.ini

This looks for ntuple files for each of the event types specified in the `event` config (you can select/deselect event types with the active/inactive fields) and produces a new set of files containing only the relevant trees for the analysis (energy, position, fit validity, pulse shape discrimination values, and true neutrino energy and reactor name for reactor neutrinos). Be warned, if using all event types, this will take some time! (of order several hours).

This app also divides the three (#alpha,n) processes into different files. They are all simulated at once, and so each of the three (#alpha,n) processes gets events from the same files. For each (#alpha,n) process, we loop over all (#alpha,n) events, and select the events of that specific (#alpha,n) process using the reconstructed energy.

<h3>LLH Scans</h3>

Right, we've now pruned our trees and can use them to build PDFs and a simulated dataset(s). In theory, we're now ready to launch a fit! But let's be steady Eddies, Cautious Carols, and Nervous Nigels, and make sure everything is behaving as intended. Running fits can be computationally and storage expensive, so it would be a shame to find out after running them that something had previously gone wrong.  

One way of doing this is with likelihood scans. The scan apps prepare the same datasets, systematics and likelihoods as the corresponding fit. They vary one parameter at a time across its configured range, holding the others at nominal values. The scan points are shifted so that the nominal value is one of them. The number of points comes from `llh_scan_points` in the fit config, or defaults to 150. Each point records the change in negative log likelihood from the nominal point. These are slices through the nominal point, not profile likelihood scans.

For the 'true' Asimov dataset, these scans should all minimise at (1,0) (where the x axis is parameter value / nominal value) i.e., by changing a parameter, you can't produce a dataset more similar to the target dataset than by using the exact values used to produce the target dataset. For the other forms of Asimov dataset, most parameters will minimise close to 1, but probably not exactly. The event types with higher rates will be closer to 1 as changes in these will have a greater impact on the LLH, and the small fluctuations causing the differences between the PDFs and Asimov dataset will be smaller proportionally.  

There are two likelihood scan apps. `fixedosc_llhscan` scans the non-oscillation parameters using the prepared likelihood, then scans each oscillation parameter by rebuilding the reactor PDFs at each point. Oscillation is fixed in those PDFs rather than varied by an OXO systematic during likelihood evaluation. This version checks the model used by `fixedosc_fit`.

To run the fixedosc LLH scan:  

> ./bin/fixedosc_llhscan cfg/fit_config.ini cfg/event_config.ini cfg/pdf_config.ini cfg/syst_config.ini

You will need to have a `deltam21` parameter, and exactly one `theta12`, `sintheta12`, or `sinsqtheta12` parameter specified in the fit config. A root file `llh_scan.root` will be saved in the output directory specified in the fit config. In the file will be a plot of LLH vs parameter value for each parameter.

The other version, `llh_scan`, uses the same configuration, dataset preparation and combined likelihood as `mcmc`. For each active parameter, including oscillation parameters, it varies that parameter with all the others held at their nominal values. Use a floating-oscillation configuration such as `cfg/mcmcDefault/` to check the MCMC model.

To run the LLH scan:  

> ./bin/llh_scan cfg/fit_config.ini cfg/event_config.ini cfg/pdf_config.ini cfg/syst_config.ini

Both LLH scan apps always write `llh_scan.root` in the configured output directory. With `save_outputs = 1`, they also write input spectra under `unscaled_pdfs`, `asimov_dists`, and, for fake data, `fakedata_dists`.

<h3>Running a Fit</h3>

Now the time has come to run a fit. When we float the oscillation parameters, we need extra dimensions for the reactors PDF to contain the true neutrino energy and the reactor core it came from. As there are 483 cores, if we have 140 bins in each energy dimension, we end up with over 9 million bins! This means marginalising down to 1D to compare to data takes a lot of time. There are few ideas to get around this, but first let's run individual fits where we don't vary the oscillation parameters (so everything remains 1D and fast), and do this at all points in a grid of oscillation parameter values. We can do this for one point with:

> ./bin/fixedosc_fit cfg/fit_config.ini cfg/event_config.ini cfg/pdf_config.ini cfg/syst_config.ini

The value of the oscillation parameters should be set (and fixed) in the `fit` config file. You will need to have a `deltam21` parameter, and exactly one `theta12`, `sintheta12`, or `sinsqtheta12` parameter. This performs a `Minuit` fit with the oscillation parameters fixed at those values. It builds the same fixed-oscillation likelihood as `fixedosc_llhscan`. Detailed output files are produced when `save_outputs` is set.

`fit_result.txt` contains the fitted parameter values, and `fit_result.root` contains the vectors of parameter names, postfit values and uncertainties, and a covariance matrix.

There will also be several subdirectories produced inside the `output_directory` set in the `fit_config`. `unscaled_pdfs` will contain each of the unscaled (normalised) PDFs without any systematics applied. `asimov_dists` will contain each of the PDFs scaled to their nominal rates with nominal systematics applied. If you're running with `fake_data = 1` in the `fit_config`, `fakedata_dists` will contain each of the PDFs scaled to their fake data value rate with systematics applied at their fake data values.

`postfit_dists` contains each component at the fitted parameter values, with systematics applied, and the sum for each dataset.

The `1dlhproj` and `2dlhproj` directories are produced by `mcmc`, not by `fixedosc_fit`.

If you want to run a full fit where everything is floated at once, there is an app to run some Markov Chain Monte Carlo. First you want to make sure everything in your fit config is looking sensible. The min and max values, and sigmas for each event type are fairly well tuned, so it's advised you leave these alone unless you know what you're doing. Certainly, for your first fit testing out the code I would leave it as is. If you have good reason to change them, you probably don't need to be reading this guide! You can run a full fit by:

> ./bin/mcmc cfg/fit_config.ini cfg/event_config.ini cfg/pdf_config.ini cfg/syst_config.ini

Set `save_outputs = 1` in the fit config to save the chain, posterior histograms, fit result and spectra. With `save_outputs = 0`, only console summaries are produced.

`Shape:ReactorOscProb` selects the reactor oscillation function for OXO's ordinary `Shape` systematic. Systematic setup prepares its 3D MC histogram before creating the Shape. Set the reactor component to `dimensions=1` and `oscillated=0`: the normal dataset builder makes an unoscillated 1D PDF, and Shape applies the oscillation weights during the fit. Put the Shape's dedicated group first in that component's group list, before detector systematics. The function uses true-energy bin centres and returns the oscillated/original count ratio for each reconstructed-energy bin.

For `Shape:ReactorOscProb`, set `param_names=deltam21,theta12`, `deltam21,sintheta12`, or `deltam21,sinsqtheta12`. The angle parameter must have the same name in the fit config. `theta12` is in degrees; the other names denote sin(theta12) and sin²(theta12). Both reactor and geoneutrino systematics convert internally to sin²(theta12). Bounds, proposal widths, priors and saved samples stay in the configured representation; a uniform prior in one representation is not uniform in the others.

In practice, run many chains in parallel on Condor with `util/submitCondor.py` (see Submitting Batch Jobs below). For example, 50 chains with the default config set in `cfg/mcmcDefault/` (oscillation parameters floated through the `Shape:ReactorOscProb` systematic):

```
python util/submitCondor.py mcmc /path/to/output/mcmcDefault_run1 \
    -r /path/to/this/repo \
    -e /path/to/your/env_file.sh \
    -f cfg/mcmcDefault/fit_config.ini \
    -i cfg/mcmcDefault/event_config.ini \
    -p cfg/mcmcDefault/pdf_config.ini \
    -s cfg/mcmcDefault/syst_config.ini \
    -n 50 \
    -w 86400 \
    -m 4000 \
    -j mcmcDefault
```

Each job writes to its own subdirectory, `mcmcDefault_run1/mcmcDefault_<i>/`, and the script sets `output_directory` in the copied fit config for you. `-w` is the wall time in seconds and `-m` the memory in MB; 1.2M steps take about 5–11 hours depending on the node.

This will load up all the pdfs, systematics, data etc. and creates the likelihood object in the same way the `llh_scan` app does, and then runs the MCMC. For each individual chain, you’ll have a number of files and subdirectories, some of which are the same as for `fixedosc_fit`.

`auto_correlations.txt` records the sample autocorrelations. `fit_result.txt` contains the result from the best step. Both are written when `save_outputs = 1`.

There will also be several subdirectories produced inside the `output_directory` set in the `fit_config`. `unscaled_pdfs` will contain each of the unscaled (normalised) PDFs without any systematics applied. `asimov_dists` will contain each of the PDFs scaled to their nominal rates with nominal systematics applied. If you're running with `fake_data = 1` in the `fit_config`, `fakedata_dists` will contain each of the PDFs scaled to their fake data value rate with systematics applied at their fake data values.

`postfit_dists` contains each component at the best sampled parameter values, with systematics applied, and the sum for each dataset.

In `1dlhproj` and `2dlhproj`, you have one- and two-parameter posterior histograms from the chain. The configured burn-in steps are excluded.

There will also be a ROOT file in the output directory containing the chain tree. Its filename is formed from the output directory and its parent. Each entry represents a step in the Markov Chain, with branches for parameter values, LLH, step time and acceptance. This is useful for MCMC diagnostics and debugging.

<h2>Submitting Batch Jobs</h2>

<h3>Individual Jobs</h3>

Every batch system will be different, but for submitting to a Condor based queue system there is a script, `util/submitCondor.py`, based on one for submitting RAT jobs originally from Josie (I think). You can use this for any exectuable. For running full MCMC fits, it’s advisable to run multiple identical fits at once in parallel as you probably want around 1 million steps.   

You can submit N jobs with:  

> python util/submitCondor.py exec_name output_dir -r /path/to/this/repo/ -e /path/to/env/file/ -f cfg/fit_config.ini -i cfg/event_config.ini -p cfg/pdf_config.ini -s cfg/syst_config.ini -n N -w walltime

The optional `-e` file must set up the complete runtime environment, including ROOT, RAT, OXO and any project-specific variables. Jobs source only this file; If `-e` is omitted, jobs use the inherited submission environment. The output directory will be created, with several subdirectories:

`error`, `output` and `log` contain the console outputs and batch log. `sh` and `submit` contain the files used to submit these jobs, and `cfg` contains all the config files used to run the fit. The `output_directory` set in the `fit_config` file will be overwritten by the command line `output_dir` argument.

This script supports the fit executables, `prune_trees` (using `-i` and `-f`), and `makeFixedOscTree` (using `-f`). Run `make_reactor_json` directly with an output JSON filename.

<h3>Fixed Oscillation Parameters Fits</h3>

First, you may want to do a grid scan of oscillation parameters, running a fit at each step. In these fits, the oscillation parameters are fixed at the values for that point in the scan, and everything else is floated in a Minuit fit. We would normally do ~500 steps for each oscillation parameter, so 250,000 individual fits. Now, you can try submitting 250,000 jobs at once but do so at your own (and your friendly neighbourhood sys-admin's) peril! Instead, we do one job for each value of one of the oscillation parameters, so that's 500 jobs each running 500 sequential fits. This number can be tuned for efficient running on your own cluster. The -d option allows you to change the number of fits ran per job.

There is a script, similar to `util/submitCondor.py` but for submitting these fixed oscillation parameter fits. This is in `util/submitFixedOscJobs.py`.

You can submit the jobs with:

> python util/submitFixedOscJobs.py fixedosc_fit output_dir -r /path/to/this/repo/ -e /path/to/env/file/ -f cfg/fit_config.ini -i cfg/event_config.ini -p cfg/pdf_config.ini -s cfg/syst_config.ini -d fits_per_job -w walltime -m requested_memory -j jobname

In your fit config, be sure to have `fake_data=1` (and `asimov=0`), and set the oscillation parameter's fake data values to be their nominal values. This way, for every step in the grid, the fake dataset will be produced with the nominal values, and this will be the same for each fit. The python script will update the nominal values of the oscillation parameters in the config file to the current values in the fixed oscillation parameters scan, so we always fit the same fake data, but with different fixed oscillation parameter values each time. 

`output_dir` will contain all the configs and logs for each fit, and the outputted root files will be produced in independent directories (one for each fit) inside that directory. The `output_directory` set in the `fit_config` file gets updated to be the one set as the command line argument.

A similar script, `util/submitFixedGeoJobs.py` runs jobs fits at 2D grid points in the U and Th geoneutrino rates. At each point a fit is run where the oscillation and geoneutrino parameters are fixed, and everything else float. Some configs exist in `cfg/geoFixStudy2026` to run the at the best fit oscillation parameters from the 2026 reactor analysis, and scan through the geoenutrino rates. The results of this can be plotted with `plotting/plotFixedGeoLLH.py`.
 
<h2>Postfit Analysis</h2>

<h3>Fixed Oscillation Fits</h3>

A lot of the postfit analysis will be done using scripts that reside in `util` (not to be confused with `src/util`), and `plotting`.

<h4>makeFixedOscTree</h4>

The first thing you'll want to do after running a set of fixed oscillation parameter fits is run `exec/makeFixedOscTree.cc`. This is currently the only postfit analysis step that requires compiling (which will be handled by the usual `Makefile`). It loops over the outputs of all the fits, and fills a tree. Each entry in the tree represents one fit, and each branch is a fit parameter. The nominal values and prefit constraints are also saved in vectors. You can run it with:

> ./bin/makeFixedOscTree cfg/fit_config.ini

The config files should be ones you've used to run one of the fits. If this takes a long time, you can submit it as a batch job using `util/submitCondor.py`:

> python util/submitCondor.py makeFixedOscTree output_dir -r /path/to/this/repo/ -e /path/to/env/file/ -f cfg/fit_config.ini -w walltime

Once the tree is made, it finds the fit that had the best LLH, and reruns that fit with the `save_outputs` bool in the fit config set to true and Minuit settings updated to do a more rigorous fit. The covariance matrix from this fit is saved along with the tree of fit results.

<h4>plotFixedOscLLH</h4>

The next thing to do is to plot the best fit LLH as a function of the oscillation parameters. You can do this by running `plotting/plotFixedOscLLH.C` over the output `TTree` from `makeFixedOscTree`. It loops over all the entries, and gets the oscillation parameter values and LLH for each fit. It then plots a 2D histogram where the X and Y axis are the oscilltion parameters, and the Z axis is the LLH. A canvas is saved in both a `.root` and `.pdf` file, in the top level output directory of the set of fits. 

Also plotted are profile LLHs for both $\Delta m^2_{21}$ and $\theta_{12}$. Each are saved as `.pdf` file and also in the same root file as the 2D plot. The $1\sigma$ bounds are calculated and saved in that root file. These can then be used by `plotFixedOscParams` as the oscillation parameter postfit uncertainties.

It can be run by doing:

> root -l 'plotting/plotFixedOscLLH.C("/path/to/makeFixedOscTree/output.root")'

<h4>plotFixedOscFullLLH</h4>

This python script makes a plot similar to `plotFixedOscLLH` but with the 1D profiles attached in the same canvas. It is run with:

> python plotting/plotFixedOscFullLLH.py -i /path/to/makeFixedOscTree/output.root

<h4>plotFixedOscDist</h4>

This script will loop over all entries in the output `TTree` from `makeFixedOscTree`, and find the fit with the minimum best LLH. It then goes to the directory of that fit, and plots the distributions saved in the `postfit_dists` directory. The distributions plotted are the data, the total MC, and each fitted component, including its systematics. Also saved is a similar plot but with components grouped together. For both plots, a panel showing the ratio of the total MC postfit prediction to the data is plotted below the main histogram.

Currently, the user hard codes the parameter names, dataset names, and axis ranges for each dataset in the top of this file. When running, the user supplies an option as the second argument: `1` plots the distributions for the first dataset, `2` plots the distributions for the second dataset, and `0` plots the summed distributions for both datasets. Obviously it would be nice if this wasn't hardcoded, and you could use N datasets, but for expedience this is how it is at the moment.

The user can also choose to supply a data directory. It will then look for files with the names `data_datasetname.root` where `datasetname` is set at the top of the script. Data from that file(s) will then be looped over and binned, with the binning also set at the top of the file. The MC distributions will then be scaled to match the 'per MeV' scaling of the data bin width.

Canvases are saved in both `.root` and `.pdf` files, in the top level output directory of the set of fits. You can run it with:

> root -l 'plotting/plotFixedOscDist.C("/path/to/makeFixedOscTree/output", dataset_option, data_directory)'

with the last two arguments being optional, defaulting to `1` and `""`.

In this script Latex labels are made for each PDF (currently reactor IBDs, Geo U, Geo Th, the three (#alpha,n)s, and Sideband). If more PDFs are used, these will be added to the backgrounds group without a Latex label. It is recommended this script gets updated if the fit parameters change. Also, the script assumes you have a TStyle called `snoplus` defined in your `rootlogon.C` file.

<h4>plotFixedOscParams</h4>

This script loops over all entries in the output `TTree` from `makeFixedOscTree`, and finds the fit with the minimum best LLH. It then plots each parameter in that entry, relative to it's nominal value. Any prefit constraints are also plotted, relative to nominal values. A canvas is saved in both a `.root` and `.pdf` file, in the top level output directory of the set of fits. The constraint, nominal, and postfit histograms are also saved in the root file. If `plotFixedOscLLH` has already been run, it will look for the outputted `LLH.root` file to find the oscillation parameter uncertainties. Otherwise, it will use the grid spacing.

The reactor parameter names are hardcoded so that they can be multiplied by the correct reactor ratio. The order the parameters are plotted in are also hardcoded.

There are two optional bool arguments the user can provide, to specify whether or not to plot correlations (defaults to true), and if the alpha-n classifier was used (defaults to false). If the former is true, a separate canvas is also saved, containing a plot of the postfit correlation matrix of all fit parameters.

You can run it with:

> root -l 'plotting/plotFixedOscParams.C("/path/to/makeFixedOscTree/output.root", plot_correlations_bool, alphan_classifier_bool)'

<h3>Prefit Plots</h3>

<h4>plotLLHScans</h4>

This script loops through the objects in a LLH scan output file and prints each histogram to a canvas. Each canvas is saved as one page in a PDF file, and in a root file. The output file is saved in the same location as the input file. You can run it with:

> root -l 'plotting/plotLLHScans.C("/path/to/llhscan/llhscan.root")'

<h4>plotLLHScanPairs</h4>

This script loops through the objects in a LLH scan output file and prints each histogram to a canvas. If a parameter exists for two datasets, the two scans for those parameters are plotted on the same canvas. If the parameter only exists for one dataset, it is plotted alone. Each canvas is saved as one page in a PDF file, and in a root file. The output file is saved in the same location as the input file. You can run it with:

> root -l 'plotting/plotLLHScanPairs.C("/path/to/llhscan/llhscan.root")'

<h4>compare2LLHScans</h4>

This script loops through the objects in the first of two LLH scan output files and prints each histogram to a canvas. If a histogram with the same name exists in the second LLH scan output files, these are drawn on the same canvas. Each canvas is saved as one page in a PDF file, and in a root file. You can run it with:

> root -l 'plotting/compare2LLHScans.C("/path/to/llhscan1/llhscan.root", "/path/to/llhscan2/llhscan.root", "fit 1 label", "fit 2 label")'

<h4>compare3LLHScans</h4>

This script is the same as `compare2LLHScans` but for three scans. It loops through the objects in the first LLH scan output file and prints each histogram to a canvas. If a histogram with the same name exists in the second and/or third LLH scan output files, these are drawn on the same canvas. Each canvas is saved as one page in a PDF file, and in a root file. You can run it with:

> root -l 'plotting/compare3LLHScans.C("/path/to/llhscan1/llhscan.root", "/path/to/llhscan2/llhscan.root", "/path/to/llhscan3/llhscan.root", "fit 1 label", "fit 2 label", "fit 3 label")'

<h4>plot1DPDFs</h4>

This script loops through all the unscaled PDF files in a directory, and plots each on a different page of a PDF file, and in a root file. You can run it with:

> root -l 'plotting/plot1DPDFs.C("/path/to/dir/")'

<h4>compare1DPDFs</h4>

This script loops through all the unscaled PDF files in a directory, and plots each on a different page of a PDF file, and in a root file. The user inputs the two directory paths, and labels for the legend. You can run it with:

> root -l 'plotting/compare1DPDFs.C("/path/to/dir1/", "/path/to/dir2/", "label1", "label2")'

<h4>compare1DDatasetPDFs</h4>

This script loops through all the unscaled PDF files in a directory, and finds the files with the siffix for one dataset, and looks for files for the same event type for the second dataset. If they both exist, it plots them both on the same page, and saves each plot on a different page of a PDF file, and in a root file. The user inputs the directory path, two dataset suffixes, and labels for the legend. You can run it with:

> root -l 'plotting/compare1DDatasetPDFs.C("/path/to/dir/", "dataset_suffix1", "dataset_suffix2", "label1", "label2")'

<h4>plot2DPDFs</h4>

This script loops through all the 2D unscaled PDF files in a directory, and plots each (with colz) on a different page of a PDF file, and in a root file. You can run it with:

> root -l 'plotting/plot2DPDFs.C("/path/to/dir/")'

<h4>projectPDFs</h4>

This script loops through all the 2D unscaled PDF files in a directory, and plots the projection of each on the y-axis (likely to be the (#alpha,n) classifier axis) on a different page of a PDF file, and in a root file. You can run it with:

> root -l 'plotting/projectPDFs.C("/path/to/dir/")'

<h3>Fit Comparisons</h3>

<h4>compareFixedOscParams</h4>

This script opens the outputted root file from `plotFixedOscParams` for two different fits, and plots the postfit parameters for each and the prefit nominals and constraints from the first. The user also inputs labels for each fit to be printed into the legend, along with the name (without file type suffix) of the outputted files. You can run it with:

> root -l 'plotting/compareFixedOscParams.C("/path/to/fit1/params.root", "fit 1 label", "/path/to/fit2/params.root", "fit 2 label", "/path/to/output/file")'

<h4>compare3FixedOscParams</h4>

This is just like `compareFixedOscParams` but it runs over 3 fits:

> root -l 'plotting/compare3FixedOscParams.C("/path/to/fit1/params.root", "fit 1 label", "/path/to/fit2/params.root", "fit 2 label", "/path/to/fit3/params.root", "fit 3 label", "/path/to/output/file")'.

<h4>compare3FixedOscDists</h4>

This script compares the total fitted distributions for one named dataset across three fits. It reads `data_<dataset>.root` from the first fit and `postfit_dists/postfitdist_<dataset>.root` from each fit, then plots their ratios to the first fit's data. It saves `distcomp_<dataset>.pdf` and `.root` in the first fit's `plots` directory.

> root -l 'plotting/compare3FixedOscDists.C("/path/to/fit1/fit_result_tree.root", "/path/to/fit2/fit_result_tree.root", "/path/to/fit3/fit_result_tree.root", "fit 1 label", "fit 2 label", "fit 3 label", "dataset1")'

<h4>compareContours</h4>

This script uses the files created by `plotFixedOscLLH` and draws the 1 $\sigma$ contours for two fits. The user inputs the two file paths and labels for the legend, along with the name (without file type suffix) of the outputted files.

> root -l 'plotting/compareContours.C("/path/to/fitdir1/plots/LLH.root", "label1", "/path/to/fitdir2/plots/LLH.root", "label2", "/path/to/output/file")'

<h4>compare3Contours</h4>

This script is just like `compareContours`, but compares three contours instead of two:

> root -l 'plotting/compare3Contours.C("/path/to/fitdir1/plots/LLH.root", "label1", "/path/to/fitdir2/plots/LLH.root", "label2", "/path/to/fitdir3/plots/LLH.root", "label3", "/path/to/output/file")'

<h4>compare3Contours2Sigmas</h4>

This script is just like `compare3Contours`, but compares three contours both at the 1 and 2 sigma levels:

> root -l 'plotting/compare3Contours2Sigmas.C("/path/to/fitdir1/plots/LLH.root", "label1", "/path/to/fitdir2/plots/LLH.root", "label2", "/path/to/fitdir3/plots/LLH.root", "label3", "/path/to/output/file")'

<h3>Postfit Utils</h3>

Inside `./util`, there are (among other things) several scripts useful for combining postfit results.

<h4>Get Integrals</h4>

There are two scripts for getting the integral of scaled PDFs over a certain range. These are particularly useful for calculating postfit event rates, where the quoted fit result may include events inside the buffer region. There are two scripts, one for looking at Asimov scaled PDFs and one comparing postfit scaled PDFs. These could obviously easily be combined into one script in the future.

They are with:

> python util/get_integrals.py /path/to/postfit_dists
> python util/get_asimov_integrals.py /path/to/asimov_dists

<h4>Combine Fit Results</h4>

This script uses the above 'get integral' scripts to get the postfit parameter event rates and proportional uncertainties, along with the postfit values and uncertainties for all other parameters, including the geoneutrino ratio. It also sums event rates over datasets to get the total number of each event rate. This is very hard coded for the parameters for the current analysis, but is useful for what it does. It outputs a CSV file of all parameter values and uncertainties. It runs `plotFixedOscParams`, so make sure `plotFixedOscLLH` has been run to get the correct oscillation parameter uncertainties. You can run it with:

> python util/combine_postfit_results.py /path/to/postfit_dists /path/to/fit_result_tree.root <correlated_fit_bool> <(alpha,n)_classifier_bool>

where the `fit_result_tree.root` was outputted by `makeFixedOscTree` and the `correlated-fit-bool` and `(alpha,n)-classifier-bool` are whether the fit used correlated normalisations and the (alpha,n) classifier. It would be nice if it automatically could determine the two bools but this is how it is for now. Also note the parameters are hard-coded inside `combine_postfit_results`. This really is meant as a handy helper for the specific current analysis than a general tool.

<h4>Run Postfit Scripts</h4>

This is a 'master script' to run nearly all of the postfit analysis steps. Because of this, there are large elements of hard-coding (including filepaths which maybe could be better read from configs), but it is very handy and worth committing. It runs `makeFixedOscTree` to combine the configured fixed oscillation fits, and then runs the plotting and fit comparison scripts. It also runs `combine_postfit_results` to make the postfit parameter value CSV file. You need to supply it with a directory that contains a `fit_config.ini` containing `reactorsjson` and the fixed-oscillation scan counts. This should be identical to one of the ones you used to run these fits, but with the `output_directory` field in the `fit_config` set to the top level output directory for the set of fits. It is run with:

> ./util/run_postfit_scripts.sh /path/to/config/dir <correlated_fit_bool> <(alpha,n)_classifier_bool> <datafit_bool>

<h4>Make Latex Tables</h4>
There are two scripts for making the style of tables in the latex report. The first has three fit results (one for each oscillation parameter prior case), along with the nominal values. The second has two fit results (one for correlated normalisations, one for uncorrelated normalisations), along with the nominal values. The inputs are CSV files, which you make by combining the outputs of `combine_postfit_results`. They are run with:

>python tables/make_prior_tables.py combined_input.csv

>python tables/make_corr_nocorr_tables.py combined_input.csv

These are very hard-coded, but at some point that's unavoidable.

The outputs are two tables for each script, one with the oscillation and normalisation parameters, and one with the systematic parameters.

There's also a script that just compares datafits with the parameter values in tables from the final paper:

>python tables/make_data_tables.py combined_input.csv


<h3>MCMC</h3>

With `save_outputs = 1`, each chain saves `fit_result.txt`, `auto_correlations.txt`, the fitted distributions, one- and two-parameter posterior histograms, and a ROOT tree of sampled steps. The configured burn-in is excluded from the posterior histograms. To combine ROOT trees from parallel chains, use `hadd`. The copied configs for batch jobs are saved by `util/submitCondor.py` in its `cfg` directory.

The older `make_plots` and `auto_corrs` executables are excluded from `make all`. Their separate Makefile targets currently fail to compile because those sources still use the old `bbfit` namespace. The MCMC output tree and the posterior histograms written by `mcmc` itself are available without these utilities.

If you've combined chain trees with `hadd`, `util/findMaxStep.C` can identify the chain with the best step:

> root -l -b -q 'util/findMaxStep.C("dataset_directory_name", "fit_name", number_of_chains)'

The `auto_correlations.txt` file written by `mcmc` contains the chain's autocorrelation values.

<h2>Physics Results</h2>

To be added!
 
