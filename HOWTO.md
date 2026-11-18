<h1>Antinu Likelihood Analysis How To Guide</h1>

This file is a guide on how to run the most commonly used apps in the antinu oxo likelihood analysis.  

The first thing you'll want to do is get the environment variable in `env.sh` pointing to your `OXO` install. It's recommended you make a copy of this `env.sh` and use that to avoid merge conflicts down the line.

All of the executables live within the `exec` directory. Within `./src`, there is code for interfacing with config files and other useful classes, methods, and functions. Templates of the config files themselves live in `./cfg`.

<h2>Source Code and Classes</h2>

<h3>Config Loaders (src/config)</h3>

For each type of config file, there is a class defined in `./src/config`, along with a loader for that class. Each of these follows a similar structure: the config class has an attribute for each field in the config file, and the config loader classes read values in from the config files and return an instance of the corresponding config class.

<h3>Configs (src/config)</h3>

There are several config files you'll need for running different apps. It's intended that you make a copy of the template config files in `./cfg` and fill in the filepaths yourself. Just be aware of any subsequent changes to the format or content of the templates, as you'll probably want to propagate these to your own local config files.

<h4>Event</h4>

This file should contain information on the types of events you want to include in the analysis. An example is shown in `cfg/event_config.ini`. There should be a section, `summary`, which controls which datasets are used in the fit, and then each dataset has it's own section which controls which event types are active and where the MC files live. Then there's a table for each event type, which contains a Latex name for axis titles, the exact location of files, the dimensions of the PDF, and the groups the PDF should belong to (different syestematics apply to different groups, see the `OXO` documentation for more info).

<h5>Summary</h5>

- `datasets`: The name of the datasets you want to fit. These should be comma separated, and each should be the name of a dataset table also defined in the file. 

<h5>Datasets</h5>

Then name of these tables should be the name of the dataset.

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
- `fit_dists`: Which PDFs should be fit. Can use `all`
- `output_directory`: Where to write output files to
- `n_steps`: Number of steps per iteration for Hamiltonian MCMC
- `epsilon`: Hamiltonian step size
- `hmc_iterations`: Number of steps in Hamiltonian MCMC
- `hmc_burn_in`: The number of steps rejected at the start of the Hamiltonian Markov Chain. This can be changed post analysis for most purposes
- `sigma_scale`: Global scaling factor each step size gets multiplied by
- `beeston_barlow`: Bool to determine if we should use the Beeston-Barlow method to account for MC stats uncertainty. If 1/True, it is used
- `save_outputs`: Bool to determine if the output files should be written. If 1/True, they are saved. This is useful for running multiple fixed oscillation fits to avoid producing ~millions of files at once. It is currently only implemented for the `fixedosc_fit` exec (for all other execs output files will be produced regardless of this bool)
- `llh_scan_points`: For the llh scan execs, this will set how many points are in the scan (defaults to 150). Other execs will just ignore this parameter
- `reactorsjson`: path to the reactor-index, name and baseline JSON. Relative paths are resolved from the executable's working directory.
- `fixedosc_dm_points`: number of fixed-fit scan points in delta-m-squared (default 500).
- `fixedosc_theta_points`: number of fixed-fit scan points in the configured theta-like parameter (default 500).

<h5>Fit Parameter Tables</h5>

The name of these tables should be the name of the fit parameter.
- `nom`: Nominal value of the parameter
- `min`: The minimum value the parameter is allowed to take
- `max`: The maximum value the parameter is allowed to take
- `sig`: The relative width of the Gaussian used to propose new step values (each gets scaled by `summary:sigma_scale`). Often referred to as the step size
- `nbins`: Number of bins used for plotting the parameter
- `tex_label`: A latex name used for the parameter
- `constraint_mean`: Mean of the prior constraint on the parameter
- `constraint_sig`: Width of the prior constraint on the parameter
- `constraint_ratiomean`: Mean of prior constraint on ratio of this parameter to another
- `constraint_ratiosigma`: Width of prior constraint on ratio of this parameter to another
- `constraint_ratioparname`: Name of other parameter the ratio constraint corresponds to. This must be another parameter defined and activated in this file.
- `constraint_corr`: Correlation of this parameter to another
- `constraint_corrparname`: Name of other parameter this parameter is correlated to. This must be another parameter defined and activated in this file. Both correlated parameters must have a `constraint_mean` and `constraint_sig` defined.
- `fake_data_val`: Value to use to produce a fake dataset if `fake_data` is true in the `Summary`
- `fix`: Whether or not parameter is fixed in the fit

A note on HMCMC: There is currently no functionality to app that runs the HMCMC. For the double beta analysis, there were studies done where some MCMC was run, followed by some HMCMC starting from the best fit points of the MCMC. However, this didn't give any improvement in results or efficiency over the straight MCMC. The functionality to read in HMCMC values is preserved in case we ever want to use it in the future, but at the moment we just use MCMC. It's likely in the future we'll just remove the functionality.

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
- `type`: The type of systematic being applied. This is not just the OXO derived systematic class, as here we also determine which function we will use for function-based systematics. The options are defined in the `SystFactory` (this will be described later in this document)
- `function`: The name of the function used if the systematic needs one. There is some redundancy here with type at the moment but this may be improved in the future
- `group`: This is the group of pdfs the systematic will be applied to. If this is left blank it will apply to all PDFs
- `dataset`: This is the dataset the systematic applies to. This should match the name of a dataset defined in the `event` config file

<h3>Data Preparation (src/data)</h3>

- `BuildDist.hh/.cc` – Builds binned PDFs from MC events using the PDF config, including flat and oscillated distributions, and constructs their axes.
- `DatasetSetup.hh` – Loads the event and PDF configs into the dataset names, active components, data paths and observables needed for a run.
- `LoadEvents.hh` – Opens pruned MC ntuples and loads observed data from a ROOT ntuple or histogram.
- `BuildDatasets.hh/.cc` – Prepares each dataset for fitting, from its MC PDFs through to the selected data distribution and initial parameter values.
- `BuildComponent.hh/.cc` – Turns prepared PDFs into expected event distributions by applying systematics and rates, then sums them for a dataset.

<h3>Fit Construction (src/fit)</h3>

- `FitParameters.hh/.cc` – Collects parameter settings from the fit config, initialises dataset values, and removes parameters and constraints unused by the selected datasets.
- `OscParams.hh` – Checks which theta12 convention is configured, converts it to sin-squared theta12, and holds the nominal oscillation values.
- `BuildLikelihood.hh/.cc` – Configures a `BinnedNLLH` for each prepared dataset, including its components, data and systematics.
- `CombineLLHs.hh` – Combines the dataset likelihoods into a `StatisticSum` and gathers their initial parameter values.
- `FitConstraints.hh` – Adds Gaussian, correlated, ratio, fractional difference and shape constraints to the combined likelihood.
- `ShapeConstraints.hh` – Defines the named functions used by shape constraints.

<h3>Systematics (src/syst)</h3>

- `SystFactory.hh/.cc` – Selects a configured systematic type and connects it to the appropriate OXO class and function.
- `SystFuncs.hh/.cc` – Implements the physics and shape functions used by systematics, including Birk's law, geoneutrino and reactor oscillation, and alpha-n classification.
- `BuildSysts.hh/.cc` – Uses `SystFactory` to create and configure one instance of each applicable systematic for every dataset, including the MC-based reactor oscillation shape.
- `ApplySysts.hh` – Applies the systematics matching a component's groups at the supplied parameter values.

<h4>Available Systematic Types</h4>

The `type` field in the syst config accepts these values (the names after a colon select a function or, for convolution, a kernel and width dependence):

- `Scale`, `Shift` – Linear axis scaling and shifting, direct from OXO.
- `Conv:Gaussian:SquareRootScale` – Gaussian smearing with a square-root-dependent width.
- `ScaleFunction:BirksLaw` – Birk's-law energy scaling
- `Shape:GeoOscProb`, `Shape:ReactorOscProb` – Geoneutrino and reactor oscillation shapes.
- `Shape:AlphaNClass{Reac,GeoU,GeoTh,BPLike,AlphaN}{PPO,BisMSB}` – Alpha-n classifier shape systematics for each listed component and scintillator.
- `Shape:reactor_nubar_norm`, `Shape:geonu_U_norm`, `Shape:geonu_Th_norm`, `Shape:alphan_PRecoil_norm`, `Shape:alphan_CScatter_norm`, `Shape:alphan_OExcited_norm`, `Shape:bipolike_norm`, `Shape:atmospheric_norm` – Shape normalisation factors.

For `Shape:ReactorOscProb`, set the reactor component to `dimensions=1` and `oscillated=0`, with the Shape's dedicated group first in its group list. `BuildSysts` uses the component's MC to prepare the oscillation shape, which is applied to the unoscillated PDF during the fit.

Set `param_names=deltam21,theta12`, `deltam21,sintheta12`, or `deltam21,sinsqtheta12`, matching the names in the fit config. `theta12` is in degrees; the other forms are sin(theta12) and sin²(theta12).

<h3>Output (src/output)</h3>

- `FitSummary.hh/.cc` – Prints a summary of the fit parameters, datasets, oscillation settings and PDF axes.
- `OutputDirs.hh` – Defines and creates the directories used for fit distributions.
- `RateReporting.hh` – Calculates component rates outside the likelihood buffer bins and prints nominal and best-fit comparisons.
- `WriteDists.hh` – Writes input PDFs, Asimov and fake-data distributions, and post-fit components and totals to ROOT files.

<h3>Utils (src/util)</h3>

- `StringUtils.hh/.cc` – Contains small string helpers for splitting delimited values and stripping quote marks.

<h2>Apps</h2>

These above classes and config files, along with all the OXO classes, are brought together in various apps in the `exec` directory. These are described in this section in the order you will probably want to use them:

<h3>Calculating Reactor Distances</h3>

To calculate the oscillation probability for an MC reactor IBD event, we need to know how far it has travelled, so we need to know the distance from the reactor it was produced in to SNO+. Rather than do this calculation on the fly at every event for every iteration of a fit, we precalculate the distance for each reactor core. You can do this by running:

> ./bin/make_reactor_json reactorjsons/allcores_rat_VERSION.json

This will loop over cores defined in the `REACTORS` RATDB table, and calculate the distance from it to SNO+. This gets saved in a JSON file, along with the reactor's name, and a unique integer assigned to the core.

The output JSON path is supplied on the command line. Set `reactorsjson` in the fit config to that path. Keep the same JSON mapping for pruning and fitting: the stored reactor indices must correspond to the MC's `reactorIndex` values. These are likely not to change very often (only when cores come online, or an error in the locations in the RATDB table is fixed. If a reactor's distance to SNO+ physically changes we probably have more pressing issues!), so they can be saved and committed in the `reacttorjsons` directory to avoid having to reproduce them. It will be good practice to indicate in the filename the RAT version it was produced with, and which cores it contains. For RAT 8.0.0, there is a committed json file with all cores, and another one with just Bruce 1 for quicker testing of the machinery.

<h3>Pruning Trees</h3>

Raw SNO+ ntuples, although much more lightweight than RATDS files, still contain much more information than we need at analysis level. So instead of carrying around that deadweight for the whole analysis, we first take the time to prune out the branches we don't need and save what we do need in new files to be used.  

The `prune_trees` app does this for us. Having set your filepaths for each dataset in your `events` config file, simply run:  

> ./bin/prune_trees cfg/event_config.ini cfg/fit_config.ini

This looks for ntuple files for each of the event types specified in the `event` config (you can select/deselect event types with the active/inactive fields) and produces a new set of files containing only the relevant trees for the analysis (energy, position, fit validity, pulse shape discrimination values, and true neutrino energy and reactor name for reactor neutrinos). Be warned, if using all event types, this will take some time! (of order several hours).

This app also divides the three $(\alpha,n)$ processes into different files. They are all simulated at once, and so each of the three $(\alpha,n)$ processes gets events from the same files. For each $(\alpha,n)$ process, we loop over all $(\alpha,n)$ events, and select the events of that specific $(\alpha,n)$ process using the reconstructed energy.

<h3>LLH Scans</h3>

Right, we've now pruned our trees and can use them to build PDFs and a simulated dataset(s). In theory, we're now ready to launch a fit! But let's be Steady Eddies, Cautious Carols, and Nervous Nigels, and make sure everything is behaving as intended. Running fits can be computationally and storage expensive, so it would be a shame to find out after running them that something had previously gone wrong.  

One way of doing this is with likelihood scans. The likelihood scan apps first build the same test statistic that the fit will use, and we will hand it the same PDFs, parameters, and systematics we intend to hand the fitter. For a single floated parameter, it will scan over a range of values centred at the nominal value. For each of 150 (or however many are set in the config) steps, it will rebuild the dataset by scaling the PDFs and applying any systematics all set at their nominal values, except the parameter in question which is set to the value we've reached in the scan. This is compared to the actual dataset using the test statistic (probably the binned LLH). The LLH (minus the nominal LLH) is saved at each step, and once we reach the end of the scan for a parameter, it is set back to its nominal value, and we repeat for the next parameter.  

For the Asimov dataset, these scans should all minimise at (1,0) (where the x axis is parameter value / nominal value) i.e., by changing a parameter, you can't produce a dataset more similar to the target dataset than by using the exact values used to produce the target dataset.  

There are two likelihood scan apps. The first, `fixedosc_llhscan` initially loops over all parameters except the oscillation parameters. The oscillation parameters are not handed to the `OXO` `BinnedNLLH` object, and the reactor PDF is produced with the oscillation parameters at nominal values. After performing the LLH scan for all non-oscillation parameters, it rebuilds the reactor PDF for each of the 150 values for each oscillation parameter. For each of these, it rebuilds the dataset and hands this to the `BinnedNLLH` to compare to the Asimov dataset. This is because, for speed, sometimes we want to run fits where the oscillation is handled outside of `OXO` (see the fixed oscillation parameters fits below). This version of the likelihood scan is designed to be used for validating that fit.

To run the fixedosc LLH scan:  

> ./bin/fixedosc_llhscan cfg/fit_config.ini cfg/event_config.ini cfg/pdf_config.ini cfg/syst_config.ini

`random_seed` in `[summary]` controls random-number generation for all fitting and scanning executables. It defaults to 4357. A positive value gives repeatable results; zero chooses a positive seed once during shared analysis setup and prints it. Each oscillated PDF construction restarts a local generator with that resolved seed, so each input event receives the same random number across scan points. MCMC uses a separate generator initialised with the same seed and advances it throughout sampling. Use the same positive `random_seed` for fixed fits and scans you want to compare. The floated-oscillation shape calculation itself does not use random event selection.

You will need to have a `deltam21` parameter, and exactly one `theta12`, `sintheta12`, or `sinsqtheta12` parameter specified in the fit config. A root file `llh_scan.root` will be saved in the output directory specified in the fit config. In the file will be a plot of LLH vs parameter value for each parameter.

The other version of the likelihood scan, `llh_scan`, hands all the parameters, including the oscillation parameters, to the `BinnedNLLH`. It scans through all parameters in the same way, setting the values in the `BinnedNLLH` and evaluating the likelihood. This version of the likelihood scan is designed to be used for validating the full fits where every parameter floats at once (see below).

To run the LLH scan:  

> ./bin/llh_scan cfg/fit_config.ini cfg/event_config.ini cfg/pdf_config.ini cfg/syst_config.ini

Both LLH scan apps will also produce several subdirectories inside the `output_directory` set in the `fit_config`. `unscaled_pdfs` will contain each of the unscaled (normalised) PDFs without any systematics applied. `asimov_dists` will contain each of the PDFs scaled to their nominal rates with nominal systematics applied. If you're running with `fake_data = 1` in the `fit_config`, `fakedata_dists` will contain each of the PDFs scaled to their fake data value rate with systematics applied at their fake data values.

<h3>Running a Fit</h3>

Now the time has come to run a fit. When we float the oscillation parameters, we need extra dimensions for the reactors PDF to contain the true neutrino energy and the reactor core it came from. As there are 483 cores, if we have 140 bins in each energy dimension, we end up with over 9 million bins! This means marginalising down to 1D to compare to data takes a lot of time. There are few ways to get around this, but first let's run individual fits where we don't vary the oscillation parameters (so everything remains 1D and fast), and do this at all points in a grid of oscillation parameter values. We can do this for one point with:

> ./bin/fixedosc_fit cfg/fit_config.ini cfg/event_config.ini cfg/pdf_config.ini cfg/syst_config.ini

The value of the oscillation parameters should be set (and fixed) in the `fit` config file. You will need to have a `deltam21` parameter, and exactly one `theta12`, `sintheta12`, or `sinsqtheta12` parameter. This performs a `Minuit` fit, with the oscillation parameters fixed at those values and the Minuit settings in the `fitCofig`. It will load up all the pdfs, systematics, data etc. and creates the likelihood object in the same way the `llh_scan` app does, and then runs a fit. A variety of output files are produced (if the `save_outputs` in the `fitConfig` is set).

`fit_results.txt` contains each parameter value for the maximum LLH point, and `fit_results.root` contains the vectors of parameter names, postfit values and uncertainties, and a covariance matrix.

There will also be several subdirectories produced inside the `output_directory` set in the `fit_config`. `unscaled_pdfs` will contain each of the unscaled (normalised) PDFs without any systematics applied. `asimov_dists` will contain each of the PDFs scaled to their nominal rates with nominal systematics applied. If you're running with `fake_data = 1` in the `fit_config`, `fakedata_dists` will contain each of the PDFs scaled to their fake data value rate with systematics applied at their fake data values.

`postfit_dists` contains the PDFs for each event type, scaled by the event type parameter value at the maximum LLH step in this chain, with each systematic’s value at the maximum LLH step applied. Also saved is the sum of these.

If you want to run a full fit where everything is floated at once, there is an app to run some Markov Chain Monte Carlo. First you want to make sure everything in your fit config is looking sensible. The min and max values, and sigmas for each event type are fairly well tuned, so it's advised you leave these alone unless you know what you're doing. Certainly, for your first fit testing out the code I would leave it as is. If you have good reason to change them, you probably don't need to be reading this guide! You can run a full fit by running:

> ./bin/mcmc cfg/fit_config.ini cfg/event_config.ini cfg/pdf_config.ini cfg/syst_config.ini

Set `save_outputs = 1` in the fit config to save the raw chain, fit result and spectra. With `save_outputs = 0`, only console summaries are produced.

In practice, run many chains in parallel on Condor with `util/submitCondor.py` (see Submitting Batch Jobs below). For example, 50 chains with the default config set in `cfg/mcmcDefault/`:

```
python util/submitCondor.py mcmc /path/to/output/mcmcDefault_run \
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

For reactor-analysis runs, use the data and Asimov configs in `cfg/reacAnalysis2025_mcmc_validation/` as templates. Proposal widths are stored directly in each free parameter's `sig`; the actual Gaussian proposal width is `sig * sigma_scale`, with `sigma_scale` set in `[summary]`. The current global scale is 0.42 and both recoil energy scales have `sig = 0.05`. To prepare another run, copy the config directory and adjust `iterations` and `output_directory` as needed; Condor submission sets the output directory in each job's copied config.

Each chain saves a ROOT `posteriors` tree containing every step, including burn-in and repeated states after rejected proposals. There is no `burn_in` setting in the fit config; select the cut afterwards with `util/mcmc_postfit.py --burn-in N`. This tree is the primary input to posterior analysis. `Step` counts proposals within the original chain; `Accepted` records whether each proposal was accepted; `LogL` is minus the NLL including constraints, so larger is better.

`random_seed` in `[summary]` is shared by fitting and scanning; the MCMC templates set it to zero to request automatic selection. If omitted it defaults to 4357. `submitCondor.py` writes a positive seed into each copied MCMC config. With an explicit nonzero seed, job `i` uses that seed plus `i` (wrapping within the positive 32-bit range). Seeds and chain identities are saved in `mcmc_metadata.json` and a `mcmc_metadata` object inside the chain file, together with ROOT and source versions, and the four configs.

MCMC saves the complete raw chain and does not retrieve or write the legacy posterior projections or likelihood autocorrelation text file. There is no projection config option. Posterior summaries, distribution plots and autocorrelations are calculated from the raw chains by the postfit tools.

The fitter still saves input and best-sample spectra: `unscaled_pdfs`, `asimov_dists`, optional `fakedata_dists`, and `postfit_dists`. For posterior analysis use the selected retained sample rebuilt by the postfit driver, rather than assuming a chain's saved best-fit spectrum corresponds to the best retained entry.

<h2>Submitting Batch Jobs</h2>

<h3>Individual Jobs</h3>

Every batch system will be different, but for submitting to a Condor based queue system there is a script, `util/submitCondor.py`, based on one for submitting RAT jobs originally. You can use this for any exectuable. For running full MCMC fits, it’s advisable to run multiple identical fits at once in parallel as you probably want around 1 million steps.   

You can submit N jobs with:  

> python utils/submitCondor.py exec_name output_dir -r /path/to/this/repo/ -e /path/to/env/file/ -f cfg/fit_config.ini -i cfg/event_config.ini  -p cfg/pdf_config.ini -s cfg/syst_config.ini  -n N -w walltime

The optional `-e` file must set up the complete runtime environment, including ROOT, RAT, OXO and any project-specific variables. Jobs source only this file; If `-e` is omitted, jobs use the inherited submission environment. The output directory will be created, with several subdirectories:

`error`, `output` and `log` contain the console outputs and batch log. `sh` and `submit` contain the files used to submit these jobs, and `cfg` contains all the config files used to run the fit. The `output_directory` set in the `fit_config` file will be overwritten by the command line `output_dir` argument.

This script supports the fit executables, `prune_trees` (using `-i` and `-f`), and `makeFixedOscTree` (using `-f`). Run `make_reactor_json` directly with an output JSON filename.

<h3>Fixed Oscillation Parameters Fits</h3>

First, you may want to do a grid scan of oscillation parameters, running a fit at each step. In these fits, the oscillation parameters are fixed at the values for that point in the scan, and everything else is floated in a Minuit fit. We would normally do ~500 steps for each oscillation parameter, so 250,000 individual fits. Now, you can try submitting 250,000 jobs at once but do so at your own (and your friendly neighbourhood sys-admin's) peril! Instead, we do one job for each value of one of the oscillation parameters, so that's 500 jobs each running 500 sequential fits. This number can be tuned for efficient running on your own cluster. The `-d` option allows you to change the number of fits ran per job.

There is a script, similar to `utils/submitCondor.py` but for submitting these fixed oscillation parameter fits. This is in `util/submitFixedOscJobs.py`. 

You can submit the jobs with:

> python utils/submitFixedOscJobs.py fixedosc_fit output_dir -r /path/to/this/repo/ -e /path/to/env/file/ -f cfg/fit_config.ini -i cfg/event_config.ini -p cfg/pdf_config.ini -s cfg/syst_config.ini -d fits_per_job -w walltime (s) -m requested_memory (MB) -j jobname

In your fit config, be sure to have `fake_data=1` (and `asimov=0`), and set the oscillation parameter's fake data values to be their nominal values. This way, for every step in the grid, the fake dataset will be produced with the nominal values, and this will be the same for each fit. The python script will update the nominal values of the oscillation parameters in the config file to the current values in the fixed oscillation parameters scan, so we always fit the same fake data, but with different fixed oscillation parameter values each time. 

`output_dir` will contain all the configs and logs for each fit, and the outputted root files will be produced in independent directories (one for each fit) inside that directory. The `output_directory` set in the `fit_config` file gets updated to be the one set as the command line argument.

A similar script, `util/submitFixedGeoJobs.py` runs jobs fits at 2D grid points in the U and Th geoneutrino rates. At each point a fit is run where the oscillation and geoneutrino parameters are fixed, and everything else float. Some configs exist in `cfg/geoFixStudy2026` to run the at the best fit oscillation parameters from the 2026 reactor analysis, and scan through the geoenutrino rates. The results of this can be plotted with `plotting/plotFixedGeoLLH.py`.
 
<h2>Postfit Analysis</h2>

Choose the command for the type of fit:

| Fit type | Command |
| --- | --- |
| Floated oscillation Minuit | `util/minuit_postfit.sh FIT_DIRECTORY` |
| MCMC | `python3 util/mcmc_postfit.py CHAIN_DIRECTORY --burn-in N` |
| Fixed-oscillation grid | `util/fixedosc_postfit.sh CONFIG_DIRECTORY CORRELATED ALPHAN DATA` |

The Minuit and MCMC paths share `plotFitParams.C` for parameters/correlations and `plotFitSpectra.C` for energy spectra. `FitParamDrawing.hh` supplies drawing code shared with the fixed-grid parameter macro. Numerical helpers for Minuit fits are in `util/minuit_results.py`; MCMC calculations are in `util/mcmc_results.py` and its plotting functions are in `util/mcmc_plots.py`.

<h3>Floated Oscillation Minuit Fits</h3>

After a floated Minuit fit has completed with `save_outputs = 1` run:

```bash
source /home/parkerw/Software/env-dev_cpp17.sh
source env_my.sh
./util/minuit_postfit.sh /path/to/output/dir
```

The fit directory must contain `fit_result.root`, `fit_result.txt`, the saved `cfg/fit_config.ini`, `cfg/event_config.ini` and `cfg/pdf_config.ini`, `postfit_dists/<component>_<dataset>.root`,
`postfit_dists/postfitdist_<dataset>.root` and `data_<dataset>.root`.
The saved configs must describe the fit that produced these files. The reader uses the supplied directory, ignoring any old `output_directory` in its config.
The standard distribution plots currently support the PPO (`2p2ppo`) and bis-MSB (`bismsb`) energy spectra and component groups. Both oscillation parameters must be free, and reactor oscillation must be applied through the shape systematic rather than baked into the PDFs.

The wrapper validates the saved fit status, parameter/error dimensions, covariance and spectra before preparing outputs. It does not rerun the fit or modify its configs, original result files or saved spectra. Rerunning the command replaces derived outputs and appends to `minuit_postfit.log`.

All products are written alongside the fit on the data disk:

| Output | Contents |
| --- | --- |
| `postfit.root` | Plotting metadata: fitted values, errors, covariance, nominals, constraints, labels and dataset/component names|
| `plots/params.pdf`, `plots/params.root` | Prefit/postfit parameter plot; ROOT histograms named `nominal`, `constraints` and `postfit` |
| `plots/correlations.pdf` | Parameter correlation plot; its histogram is also saved as `correlations` in `params.root` |
| `plots/dist_all.pdf`, `plots/distgroup_all.pdf` | Combined component and grouped energy spectra, when both datasets are present |
| `plots/dist_2p2ppo.pdf`, `plots/distgroup_2p2ppo.pdf` | PPO spectra |
| `plots/dist_bismsb.pdf`, `plots/distgroup_bismsb.pdf` | bis-MSB spectra |
| Corresponding distribution `.root` files | ROOT canvases for the energy plots |
| `postfit_summary.csv` | Parameters, component event counts and U/Th ratios in a table |
| `postfit_central.json`, `postfit_uncert.json` | Numerical summary at full saved precision |
| `postfit_summary_metadata.json` | Fit source, NLLH, integration-bin convention and uncertainty method |
| `minuit_postfit.log` | Output from preparation and plotting |

The numerical helpers for floated Minuit and fixed-oscillation fits are in `util/minuit_results.py`. It prepares metadata and calculates the summary independently of plotting. `plotting/plotFitParams.C` and the fixed oscillation grid macro `plotting/plotFixedOscParams.C` share their parameter and
correlation drawing code in `plotting/FitParamDrawing.hh`. The distribution macro `plotting/plotFitSpectra.C` accepts the fit directory directly for this workflow. Observed data are rebinned from the saved data spectrum onto the 0.4 MeV display bins, whose edges must align with the saved bins. The MC is scaled to events per display-bin width. Asimov/fake-data plots retain the fine binning and MC/data ratio panel. Raw data ntuples are not reopened.

To redraw plots after running `minuit_postfit.sh`, use:

```bash
root -l -b -q 'plotting/plotFitParams.C("/path/to/fit/postfit.root")'
root -l -b -q 'plotting/plotFitSpectra.C("/path/to/fit/postfit.root",0,"","/path/to/fit",true)'
```

The plotting commands reuse `postfit.root`, which `minuit_postfit.sh` prepares automatically in the fit directory. `plotFitParams.C` is the shared parameter and correlation plotting macro for Minuit and MCMC. Its second argument defaults to `false` for Minuit errors; MCMC passes `true` for asymmetric posterior intervals. The last argument in the distribution command is `true` for data and `false` for Asimov/fake data. The dataset choice is `0` for combined, `1` for PPO and `2` for bis-MSB. The wrapper uses the fit directory for all derived outputs.

Parameter errors and correlations use the joint Minuit Hessian covariance. Boundary parameters and multimodal oscillation likelihoods may need profile studies; Hessian errors do not describe disconnected allowed regions.
Event counts are integrated over the energy likelihood bins specified by `llh_buffer_bin_low` and `llh_buffer_bin_high` (defaults: ROOT bin 8 through `number_of_bins - 20`, which is bin 150 for the standard 170-bin spectrum).

Yield errors propagate the covariance of normalisation parameters, including shared factors between datasets and the U/Th covariance. Explicit yield changes from shape or oscillation variations are not included. These are not full uncertainties on the integrated prediction. The CSV retains the existing `Central`, `Uncertainty` and `Parameter` columns for the table tools.

The existing comparison tools accept the new `plots/params.root` directly:

```bash
root -l -b -q 'plotting/compareFitParams.C("/path/to/fixed/plots/params.root","Fixed","/path/to/floated/plots/params.root","Floated","/path/to/comparison")'
```

For three fits use `compare3FitParams.C` with three filename/label pairs and an output basename. Choose a basename and create its parent directory first. MCMC inputs use the saved `posterior_intervals` graph (median and asymmetric 16th–84th percentile errors); fixed-oscillation and Minuit inputs retain their saved fitted values and errors. Comparisons check that parameter labels and ordering match, and expand the vertical range to include the intervals.

Tables are a separate comparison step, for example:

```bash
python3 tables/make_old_new_fit_tables.py /path/to/fit1/postfit_summary.csv /path/to/fit2/postfit_summary.csv \
    --fit3 /path/to/fit3/postfit_summary.csv \
    --label1 "Label 1" --label2 "Label 2 fixed" \
    --label3 "Floated" --name "data alphan" --totals-only --table-type osc \
    -o /path/to/comparison.tex
```

Use `--no-nominal` for a table containing only the parameter name and the compared fits. By default the nominal column comes from the first fit: parameter `nom` values in its saved config and physical counts from its nominal spectra, integrated outside the configured buffer bins. MCMC postfit outputs supply nominal counts through `postfit.root`. The comparison script accepts MCMC output directories or their `postfit.root` files directly; existing CSV inputs remain supported. These values are shown without an uncertainty; proposal widths and postfit errors are not nominal uncertainties. Use `--nominal-dir /path/to/fit` to select another fit as the nominal source, or `--nominal-csv /path/to/table.csv` to explicitly use a historical nominal table and its quoted uncertainties. CSV-only comparisons need a nominal source or `--no-nominal`. No historical analysis CSV is selected automatically; the caption `--name` never selects numerical inputs. Omit `--totals-only` to retain dataset-specific rows and omit `--table-type osc` to include the systematic table too.

A single Minuit result does not provide a likelihood surface or profile curves. The wrapper therefore does not run `makeFixedOscTree`, `plotFixedOscLLH.C` or `plotFixedOscFullLLH.py`.

<h3>Fixed Oscillation Fits</h3>

A lot of the postfit analysis will be done using scripts that reside in `util` (not to be confused with `src/util`), and `plotting`.

<h4>makeFixedOscTree</h4>

The first thing you'll want to do after running a set of fixed oscillation parameter fits is run `exec/makeFixedOscTree.cc`. This is currently the only postfit analysis step that requires compiling (which will be handled by the usual `Makefile`). It loops over the outputs of all the fits, and fills a tree. Each entry in the tree represents one fit, and each branch is a fit parameter. The nominal values and prefit constraints are also saved in vectors. You can run it with:

> ./bin/makeFixedOscTree cfg/fit_config.ini

The config files should be ones you've used to run one of the fits. If this takes a long time, you can submit it as a batch job using `utils/submitCondor.py`:

> python utils/submitCondor.py makeFixedOscTree output_dir -r /path/to/this/repo/ -e /path/to/env/file/ -f cfg/fit_config.ini -w walltime

Once the tree is made, it finds the fit that had the best LLH, and reruns that fit with the `save_outputs` bool in the fit config set to true and Minuit settings updated to do a more rigorous fit. The covariance matrix from this fit is saved along with the tree of fit results.

<h4>plotFixedOscLLH</h4>

The next thing to do is to plot the best fit LLH as a function of the oscillation parameters. You can do this by running `plotting/plotFixedOscLLH.C` over the output `TTree` from `makeFixedOscTree`. It loops over all the entries, and gets the oscillation parameter values and LLH for each fit. It then plots a 2D histogram where the X and Y axis are the oscilltion parameters, and the Z axis is the LLH. A canvas is saved in both a `.root` and `.pdf` file, in the top level output directory of the set of fits. 

Also plotted are profile LLHs for both $\Delta m^2_{21}$ and $\theta_{12}$. Each are saved as `.pdf` file and also in the same root file as the 2D plot. The $1\sigma$ bounds are calculated and saved in that root file. These can then be used by `plotFixedOscParams` as the oscillation parameter postfit uncertainties.

It can be run by doing:

> root -l 'plotting/plotFixedOscLLH.C("/path/to/makeFixedOscTree/output.root")'

<h4>plotFixedOscFullLLH</h4>

This python script makes a plot similar to `plotFixedOscLLH` but with the 1D profiles attached in the same canvas. It is run with:

> python plotting/plotFixedOscFullLLH.py -i /path/to/makeFixedOscTree/output.root

<h4>plotFitSpectra</h4>

This is the shared spectrum plotter for fixed-oscillation grids, floated Minuit fits and reconstructed MCMC fits. For a single fit, pass its `postfit.root`, dataset choice, an empty data-directory argument, the fit directory, and whether it is observed data, as shown in the Minuit commands above. For MCMC the fit directory is the postfit output’s `dists/` directory. The grid workflow is described below.

This script will loop over all entries in the output `TTree` from `makeFixedOscTree`, and find the fit with the minimum best LLH. It then goes to the directory of that fit, and plots the distributions saved in the `postfit_dists` directory. The distributions plotted are the data, the total MC (sum of all scaled PDFs with systematics applied), and each individual PDF scaled (without systematics applied). Also saved is a similar plot but with PDFs grouped together. For both plots, a panel showing the ratio of the total MC postfit prediction to the data is plotted below the main histogram.

Currently, the user hard codes the parameter names, dataset names, and axis ranges for each dataset in the top of this file. When running, the user supplies an option as the second argument: `1` plots the distributions for the first dataset, `2` plots the distributions for the second dataset, and `0` plots the summed distributions for both datasets. Obviously it would be nice if this wasn't hardcoded, and you could use N datasets, but for expedience this is how it is at the moment.

The user can also choose to supply a data directory. It will then look for files with the names `data_datasetname.root` where `datasetname` is set at the top of the script. Data from that file(s) will then be looped over and binned, with the binning also set at the top of the file. The MC distributions will then be scaled to match the 'per MeV' scaling of the data bin width.

Canvases are saved in both `.root` and `.pdf` files, in the top level output directory of the set of fits. You can run it with:

> root -l 'plotting/plotFitSpectra.C("/path/to/makeFixedOscTree/output", dataset_option, data_directory)'

with the last two arguments being optional, defaulting to `1` and `""`.

In this script Latex labels are made for each PDF (currently reactor IBDs, Geo U, Geo Th, the three $(\alpha,n)$, and Sideband). If more PDFs are used, these will be added to the backgrounds group without a Latex label. It is recommended this script gets updated if the fit parameters change. Also, the script assumes you have a TStyle called `snoplus` defined in your `rootlogon.C` file.

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

<h4>compare1DDatasetPDFs</h4>

This script loops through all the unscaled PDF files in a directory, and finds the files with the siffix for one dataset, and looks for files for the same event type for the second dataset. If they both exist, it plots them both on the same page, and saves each plot on a different page of a PDF file, and in a root file. The user inputs the directory path, two dataset suffixes, and labels for the legend. You can run it with:

> root -l 'plotting/compare1DDatasetPDFs("/path/to/dir/", "dataset_syffix1", "dataset_syffix2", "label1", "label2")'

<h4>plot2DPDFs</h4>

This script loops through all the 2D unscaled PDF files in a directory, and plots each (with colz) on a different page of a PDF file, and in a root file. You can run it with:

> root -l 'plotting/plot2DPDFs("/path/to/dir/")'

<h4>projectPDFs</h4>

This script loops through all the 2D unscaled PDF files in a directory, and plots the projection of each on the y-axis (likely to be the (#alpha,n) classifier axis) on a different page of a PDF file, and in a root file. You can run it with:

> root -l 'plotting/projectPDFs("/path/to/dir/")'

<h3>Fit Comparisons</h3>

<h4>compareFitParams</h4>

This script opens `plots/params.root` from `plotFitParams` or `plotFixedOscParams` for two fits (fixed oscillation, floated Minuit or MCMC), and plots the postfit parameters for each and the prefit nominals and constraints from the first. The user also inputs labels for each fit to be printed into the legend, along with the name (without file type suffix) of the outputted files. You can run it with:

> root -l 'plotting/compareFitParams.C("/path/to/fit1/params.root", "fit 1 label", "/path/to/fit2/params.root", "fit 2 label", "/path/to/output/file")'

<h4>compare3FitParams</h4>

This is just like `compareFitParams` but it runs over 3 fits:

> root -l 'plotting/compare3FitParams.C("/path/to/fit1/params.root", "fit 1 label", "/path/to/fit2/params.root", "fit 2 label", "/path/to/fit3/params.root", "fit 3 label", "/path/to/output/file")'.

<h4>compare3FixedOscDists</h4>

This script is similar to `plotFitSpectra`, but plots the total postfit scaled MC for 3 fits (without each scaled component PDF). It loops through the entries in the output `TTrees` from `makeFixedOscTree`, and finds the fit with the minimum best LLH for each. The distributions plotted are the data for the first fit, and the total MC (sum of all scaled PDFs with systematics applied) for each of the three fits, and a panel showing the ratio of each total MC postfit prediction to the data for the first fit.

> root -l 'plotting/compare3FixedOscDists.C("/path/to/fit1/params.root", "/path/to/fit2/params.root", "/path/to/fit3/params.root", "fit 1 label", "fit 2 label", "fit 3 label")'

<h4>compareContours</h4>

This script uses the files created by `plotFixedOscLLH` and draws the 1 $\sigma$ contours for two fits. The user inputs the two file paths and labels for the legend, along with the name (without file type suffix) of the outputted files.

> root -l 'plotting/compareContours("/path/to/fitdir1/plots/LLH.root", "label1", "/path/to/fitdir2/plots/LLH.root", "label2", "/path/to/output/file")'

<h4>compare3Contours</h4>

This script is just like `compareContours`, but compares three contours instead of two:

> root -l 'plotting/compare3Contours("/path/to/fitdir1/plots/LLH.root", "label1", "/path/to/fitdir2/plots/LLH.root", "label2", "/path/to/fitdir3/plots/LLH.root", "label3", "/path/to/output/file")'

<h4>compare3Contours2Sigmas</h4>

This script is just like `compare3Contours`, but compares three contours both at the 1 and 2 sigma levels:

> root -l 'plotting/compare3Contours2Sigmas("/path/to/fitdir1/plots/LLH.root", "label1", "/path/to/fitdir2/plots/LLH.root", "label2", "/path/to/fitdir3/plots/LLH.root", "label3", "/path/to/output/file")'

<h3>Postfit Utils</h3>

Inside `./util`, there are (among other things) several scripts useful for combining postfit results.

<h4>Get Integrals</h4>

There are two scripts for getting the integral of scaled PDFs over a certain range. These are particularly useful for calculating postfit event rates, where the quoted fit result may include events inside the buffer region. There are two scripts, one for looking at Asimov scaled PDFs and one comparing postfit scaled PDFs. These could obviously easily be combined into one script in the future.

They are with:

> python get_integrals.py /path/to/postfit_dists
> python get_asimov_integrals.py /path/to/asimov_dists

<h4>Fixed-oscillation fit summary</h4>

This script uses the above 'get integral' scripts to get the postfit parameter event rates and proportional uncertainties, along with the postfit values and uncertainties for all other parameters, including the geoneutrino ratio. It also sums event rates over datasets to get the total number of each event rate. This is very hard coded for the parameters for the current analysis, but is useful for what it does. It outputs a CSV file of all parameter values and uncertainties. It runs `plotFixedOscParams`, so make sure `plotFixedOscLLH` has been run to get the correct oscillation parameter uncertainties. You can run it with:

> python util/fixedosc_summary.py /path/to/postfit_dists /path/to/fit_result_tree.root <correlated_fit_bool> <(alpha,n)_classifier_bool>

where the `fit_result_tree.root` was outputted by `makeFixedOscTree` and the `correlated-fit-bool` and `(alpha,n)-classifier-bool` are whether the fit used correlated normalisations and the (alpha,n) classifier. It would be nice if it automatically could determine the two bools but this is how it is for now. Also note the parameters are hard-coded inside `fixedosc_summary`. This really is meant as a handy helper for the specific current analysis than a general tool.

<h4>Run Postfit Scripts</h4>

This is a 'master script' to run nearly all of the postfit analysis steps. Because of this, there are large elements of hard-coding (including filepaths which maybe could be better read from configs), but it is very handy and worth committing. It runs `makeFixedOscTree` to combine all the 500 fixed oscillation fits, and then runs `plotFixedOscLLH`, `plotFixedOscParams`, `plotFixedOscLLH` and `plotFixedOscFullLLH` with the appropriate options to make all postfit plots. It then also runs `fixedosc_summary` to make the post-fit parameter value CSV file. You need to supply it with a directory that contains a `fit_config.ini` containing `reactorsjson` and the fixed-oscillation scan counts. This should be identical to one of the ones you used to run these fits, but with the `output_directory` field in the `fit_config` set to the top level output directory for the set of fits. It is run with:

> ./util/fixedosc_postfit.sh /path/to/config/dir <correlated_fit_bool> <(alpha,n)_classifier_bool> <datafit_bool>

<h4>Make Latex Tables</h4>
There are two scripts for making the style of tables in the latex report. The first has three fit results (one for each oscillation parameter prior case), along with the nominal values. The second has two fit results (one for correlated normalisations, one for uncorrelated normalisations), along with the nominal values. The inputs are CSV files, which you make by combining the outputs of `fixedosc_summary`. They are run with:

>python make_prior_table.py combined_input.csv

>python make_corr_nocorr_table.py combined_input.csv

These are very hard-coded, but at some point that's unavoidable.

The outputs are two tables for each script, one with the oscillation and normalisation parameters, and one with the systematic parameters.

There's also a script that just compares datafits with the parameter values in tables from the final paper:

>python make_data_table.py combined_input.csv


<h3>MCMC</h3>

Build the tools and run one driver:

```bash
source /path/to/env/file
source env_my.sh
make OXSX_LIB_DIR=/path/to/oxo/cmake-build/src all

python3 util/mcmc_postfit.py \
    /path/to/output/chains/directory \
    --burn-in 2000
```

`--output-dir` is optional: by default results go in `<input directory>/output/`. For a ROOT-file input, they go in `output/` beside that file. With multiple inputs, the first input determines the default location. Use `--output-dir /path/to/output/` to override it.

Choose the burn-in from the traces and between-chain comparisons; 2000 above is an example for short tuning chains. `Step >= burn_in` is applied separately in each chain. Accepted and rejected steps both contribute with their original multiplicity. The driver checks step continuity, expected entry counts, successful job logs when present, free branches, and likelihood config compatibility. Duplicate inputs and merged trees with resetting steps are rejected. Multiple input directories or explicit original ROOT files may be supplied.

The output includes:

- `report.txt`: a quick run summary: chains and retained steps used, the maximum-LLH step selected, convergence warnings and reconstruction timing. R-hat above 1.01 flags disagreement between chains; low effective sample size means the correlated steps contain little independent information
- `logs/`: reconstruction, autocorrelation and ROOT plotting logs (`reconstruct.log`, `autocorrelations.log`, `root_plots.log`)
- `dists/`: reconstructed nominal and best-fit spectra used by the existing distribution plotter
- `postfit.root`: parameter summaries and diagnostics, covariance and correlation matrices, best-step values and source location, nominal and reconstructed yields, and analysis settings. Optional transformed posterior draws and numerical plot histograms are stored here too
- `parameter_table.csv`: one numerical table of parameter names, medians, asymmetric uncertainties (`error_low`, `error_high`) and 16th/84th percentile bounds (`interval_low`, `interval_high`). Includes all fitted parameters except `LogL`, plus physical rates and ratios when their posterior uncertainties were calculated. Raw fitted normalisations include buffer regions; physical-rate rows use non-buffer counts. No LaTeX tables or table PDF are generated; bespoke comparison tables can read `postfit.root` directly
- `plots/posteriors.pdf`, `traces.pdf`, `autocorrelations.pdf`, `acceptance.pdf`: parameter distributions and per-chain diagnostics; traces include a reduced drawing of the pre-cut samples. Acceptance is shown separately for each chain in 40 approximately equal blocks after burn-in (or one block per step if fewer than 40 steps remain), at the mean step of each block
- `plots/params.pdf`, `params.root`, `correlations.pdf`: the shared ROOT parameter and correlation plots; the default MCMC result is the median with its asymmetric equal-tailed 68% interval. This single-fit plot also shows overlapping translucent bands for the mean with standard deviation, Gaussian mean with fitted width, and mode with its highest-density 68% region. The compatibility histogram still stores medians with posterior standard deviations, while `posterior_intervals` stores the default asymmetric result
- `plots/oscillation_posterior.pdf`: binned joint posterior and 68%/95% credible-density contours, distinct from fixed-grid likelihood confidence contours; binning detail follows the effective sample size
- With `--rate-samples N`, `plots/rate_posteriors.pdf`: physical event counts outside the likelihood buffer bins, and U/Th ratios, from the selected draws
- `plots/dist*` and `distgroup*`: the existing spectrum plots for PPO, bis-MSB and their sum, rebuilt at the selected best retained sample
- `cfg/`: copies of the input configs; input paths, analysis versions and settings are recorded inside `postfit.root`

The default MCMC extraction is the **median with the equal-tailed 68% interval**, bounded by the 16th and 84th percentiles. Parameter plots, tables and fit comparisons use this same definition. These bounds describe the sampled parameter distribution, rather than sampling error on the median itself.

`posteriors.pdf` shows one pooled histogram per parameter using every retained step across all chains, including repeated states after rejected proposals. Only the median interval is shaded, with a matching shaded legend key. Mean, Gaussian fit, mode and maximum-LLH estimates remain visible as lines, with the corresponding widths described in the legend. Titles and axis labels use the saved config's `tex_label` fields, with `Steps` on the vertical axis. Legend values use two decimal places, with scientific notation for small values.

`rate_posteriors.pdf` uses the same presentation for the selected non-buffer event-count draws and U/Th ratios. Its maximum-LLH marker comes from reconstructing the actual best retained step, even when that step was not in the selected subset. Mean/median summaries and credible intervals use the original sample values. The Gaussian fit and mode use 80-bin histograms; a failed Gaussian fit is labelled unavailable.

Alternative estimates remain saved for inspection: Gaussian fitted mean with fitted width, mean with standard deviation, and mode with its highest-density 68% region. Highest-density bins are selected until they contain at least 68% of samples, with actual coverage in `hpd_mass` and possibly disconnected bounds in the `highest_density_regions` tree in `postfit.root`. These alternative bands are shown on the single-fit `params.pdf`; posterior histograms shade only the median interval. The mode legend gives its value and asymmetric distances to the highest-density bounds; for disconnected regions, these are the bounds of the region containing the mode, and the legend reports the number of regions.

The driver calls `plotting/plotFitParams.C` directly for the parameter and correlation plots. To redraw just those plots from an existing MCMC result:

```bash
root -l -b -q 'plotting/plotFitParams.C("/path/to/chains/output/postfit.root",true)'
```

It also runs `bin/auto_corrs` separately on each original chain, using the requested `--burn-in` and `--max-lag`. Results are combined in `plots/autocorrelations.pdf`, with the histograms stored inside `postfit.root`; per-chain exchange files are temporary. With `--no-plots`, `auto_corrs` is skipped, while the Python R-hat and effective sample size calculations still run.

`mcmc_postfit.py` also reports mcmc diagnostic information. It reports R-hat above 1.01 and low effective sample sizes, indicating chians haven't yet converged. Diagnostics split each chain and use equal-length prefixes if retained chain lengths differ. Posterior summaries use all retained draws. The diagnostics follow the [Stan posterior-analysis definitions](https://mc-stan.org/docs/reference-manual/analysis.html), using rank normalisation, folded split R-hat and Geyer's positive monotone sequence for ESS. Reported ESS is conservatively capped at the number of draws. Quantiles and moments use the actual samples, not smoothed or Gaussian-fitted histograms.

<h4>Rates outside the buffer regions</h4>

`util/mcmc_postfit.py` directly selects reconstruction points, invokes `bin/reconstruct_mcmc_fit`, summarises yields and writes tables. The C++ helper reconstructs spectra and physical yields, and shares `BuildComponents` and `GetComponentRates` with the fitter. It applies all sampled shapes and oscillation parameters before integrating outside the configured buffer bins. The output quantities are expected event counts for the configured exposure. Per-dataset counts, totals, and U/Th ratios are calculated sample by sample before their posterior summaries.

By default, only the highest-score retained step is reconstructed. This produces the best-fit spectra under `dists/` and physical component yields, totals and U/Th ratios in the `best_fit` tree in `postfit.root`. These best-step yields have no posterior uncertainty attached. Parameter posterior summaries and tables always use all retained draws; physical-yield rows enter the posterior tables only when yield uncertainties are explicitly requested.

Add `--rate-samples N` to calculate physical-yield posterior uncertainties using a stratified random subset of N retained draws, with selection seed 1729. Repeated states retain their correct multiplicities; duplicate parameter points are evaluated once and expanded back into the selected draws. Use `--rate-samples 0` to evaluate every retained draw. Check that intervals are stable as N increases: additional evaluations reduce subsampling error but cannot repair poor chain mixing.

`postfit.root` contains `parameter_summary` (medians, intervals, moments and diagnostics), `nominal` and `best_fit` trees, `covMatrix` and `corrMatrix`, and an `analysis` text object recording settings, best-step location, per-chain results and reconstruction timings. The optional `derived_samples` tree contains physical yields and ratios alongside the pooled sample index, chain index and original step. Posterior, autocorrelation and joint-density histograms are added to directories inside the same ROOT file when plotting. ROOT parameter and spectrum renderers are also stored under `plots/`.

Useful options:

```bash
# Fast parameter-only analysis and diagnostics
python3 util/mcmc_postfit.py BATCH --burn-in 10000 --output-dir OUTPUT --skip-rates

# Numerical summaries and tables without plotting
python3 util/mcmc_postfit.py BATCH --burn-in 10000 --output-dir OUTPUT --no-plots

# Optional physical-yield posterior uncertainties
python3 util/mcmc_postfit.py BATCH --burn-in 10000 --output-dir OUTPUT --rate-samples 2000

# Additional joint posterior plots, without repeated ROOT tree scans
python3 util/mcmc_postfit.py BATCH --burn-in 10000 --output-dir OUTPUT \
    --pairs geonu_U_norm:geonu_Th_norm p_recoil_energy_scale:class_a_ppo
```

`--pairs all` requests every free-parameter pair and can be expensive. `--max-lag` controls the diagnostic plot range, defaulting to 10000. Fixed parameters are excluded from posterior summaries and correlations.

<h4>Standalone autocorrelations</h4>

```bash
./bin/auto_corrs original_chain.root --burn-in 10000 --max-lag 10000
```

This writes `original_chain_autoCorr.root` with per-parameter histograms in `Auto_Corrs`. Use `--output path.root` to choose another destination. Lag zero is the first visible bin and equals one. The maximum lag is inclusive and is capped at the number of retained entries minus one. Constant parameters are skipped. With no options, all steps are retained and lags zero through 1000 are calculated. The zero-padded FFT uses the full zero-lag sum of squared deviations for normalisation, matching the original direct calculation. Run it before combining chains: autocorrelations must not cross chain boundaries.
