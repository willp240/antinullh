<h1>Antinu Likelihood Analysis How To Guide</h1>

This file is a guide on how to run the most commonly used apps in the double beta likelihood analysis.  

Most of the apps and source code live within the `fit` directory. Within `fit/src`, there is code for interfacing with config files. The config files themselves live in `cuts`, `pdfs`, `rates`, `results`, and `systs`.  

The first thing you'll need to do is update the path where the raw ntuples you're going to use are saved. In the `rates` config file, change `orig_base_dir` to point to the top-level directory containing them. At this point, you may also want to update `pruned_ntup_dir` to point to where you want your pruned ntuples to be outputted, `split_ntuple_fake` to point to where to output trees of pruned events used to create fake datasets, and `split_ntuple_pdf` to point to where to output trees of pruned events used to create PDFs.  

You'll also want to get the environment variables in `env.sh` pointing to your OXO install, your `bb-likelihood-analysis/fit` directory, and the top level where data will be written to. From now on, paths in this document are relative to `bb-likelihood-analysis/fit` (except for postfit outputted files).  

<h2>Pruning Trees</h2>

Raw SNO+ ntuples, although much more lightweight than RATDS files, still contain much more information than we need at analysis level. So instead of carrying around that deadweight for the whole analysis, we first take the time to prune out the branches we don't need and save what we do need in new files to be used.  

The `make_trees` app does this for us. Having set your filepaths in your `rates` config file, simply run (from within `bb-likelihood-analysis/fit`):  

>> ./bin/make_trees ../rates/config_file.cfg

This looks for ntuple files for each of the event types specified in the `rates` config (you can select/deselect event types with the active/inactive fields) and produces a new set of files containing only the relevant trees for the analysis (energy, position, fit validity, effective radius, classifiers and pulse shape descrimination values). Be warned, if using all event types, this will take some time! (Of order several hours).  

<h2>Splitting Events</h2>

Now we have our pruned events, we're going to want to make some fake datasets and PDFs. However, we may wish to make each of these from different events, to keep them statistically independent. Currently there's two ways of doing this. The simplest is to use `split_half`:  

>> ./bin/split_half ../rates/config_file.cfg fraction

This simply loops over all events of each event type and puts the first *n* in a directory for fake data, and the rest in a directory for PDFs. The paths to those directories are what you set `split_ntuple_fake` and `split_ntuple_pdf` in your `rates` config. You determine *n* by setting `fraction`, e.g., if you want half in each set `fraction` as 0.5. If you want a 1/4 in `split_ntuple_fake` and 3/4 in `split_ntuple_pdf`, use 0.25.  

Alternatively, you could use `split_data`. This goes a step further and combines events into datasets. The livetime and number of datasets produced can be set by the user:  

>> ./bin/split_data ../rates/config_file.cfg livetime(yrs) nDataSets replace useSequenceFlags

`livetime` should be the livetime of the datasets in years. `nDatasets` is the number of datasets you want to produce. `replace` is a bool. If set to `true`, once put into a dataset, the event is 'replaced' in the pool of events to be selected, so can be selected for another dataset too. If `false`, once selected, the event can't be used in a different dataset so the datasets are completely independent. Any events leftover are then saved in the `split_ntuple_pdf` directory. `useSequenceFlags` is another bool. If `true`, when selecting events, the `split_method` field is read for each event type in the `rates` config file. If this is `sequential`, events are just selected in the order they are saved in the ntuples. If it is `random`, they are selected randomly. This takes longer, but is preferable for some event types if there are correlations between consecutive events e.g., for solar events, the rate is correlated with the position of the Sun, and as the production events are saved in MC time order, selecting consecutive events will introduce a bias in the dataset.  

Running this for all event types will take < 1 hour, but longer for more datasets.  

<h2>Making the PDFs</h2>

Unless we ran `split_data` with `replace=1`, we now have a file of events for each event type in our `split_ntuple_pdf` directory. It's now time to build a PDF for each:  

>> ./bin/make_pdfs ../rates/config_file.cfg ../pdfs/config_file.cfg ../cuts/config_file.cfg

This loops over the files for each event type, produces the event distributions in the desired dimensions, having applied cuts, and then normalises each one. Be sure to have set `pdf_dir` in the PDF config to be where you've outputted the events earmarked for going into a PDF (probably the `split_ntuple_pdf` in the rates config). The cuts config contains details of the OXO analysis cuts you want to use. 

This should take < 10 minutes to run.  

<h2>Building the Asimov Dataset</h2>

For the time being, we don't have any juicy tellurium loaded data. But even if we did, we'd still want to produce simulated datasets to test out our analysis framework and determine expected sensitivities before opening the box on data. The most common form of simulated dataset is known as the Asimov dataset. This is produced by having every event type rate and systematic set to their nominal rate/value. Every floated parameter of the ultimate fit is set to the prefit expected value.  

The 'purest' Asimov dataset is built from the PDFs themselves, which are just scaled to the expected rates. This way, the dataset is exactly the sum of the scaled PDFs.  

Alternatively, you can build the dataset from the events saved in `split_ntuple_fake`. These are events from the same sample of MC but separated from those used to build the PDF so are statistically independent. This dataset will be very close to, but not exactly identical to the sum of scaled PDFs. This dataset would reflect the data we would expect if nature was exactly as we predict (all rates at nominal values) but measured non-integer numbers of events.  

Finally, you can build the dataset from events separate to the PDF events (as above), but only have integer number of events in each bin. This is done by taking a Poisson fluctuation around the expected numbers. This reflects what we expect to see in the detector if our model is identical to nature. In general, if doing this we would produce multiple datasets, to cover the range of differences we'll expect to see from statistical fluctuations. These datasets are what are produced when running `split_data`.  

To build either of the first two types of Asimov dataset, run `build_azimov`:  

>> ./bin/build_azimov ../rates/config_file.cfg ../pdfs/config_file.cfg ../cuts/config_file.cfg ../systs/config_file.cfg livetime output_file loadPDF generated_scale loading_scale

The `systs` config file contains details of the systematics you want to apply to the dataset. `livetime` is the livetime in years. `output_file` is where you want to save the dataset. This should have no file extension, a `.root` and `.h5` file will be produced. `loadPDF` is a bool, if `true`, the 'true' Asimov dataset described above will be produced. If `false`, the second option above will be produced. `generated_scale` is the amount you have scaled the the number of generated events by to get the sample of events you are building this dataset from e.g., if you ran `split_half` with `fraction=0.5`, this should be 0.5. `loading_scale` is how much you want to scale up the event rates of certain event types to approximately replicate the effect of increased tellurium loading. The event types you want to scale up with this factor are set in the `rates` config file by having `scales_with_loading=true`.  

To run `build_azimov`, first you will have had to run `split_half` rather than `split_data` as `split_data` would have already combined all the `split_ntuple_fake` into a fakedaset.  

This should take < 10 minutes to run. 

<h2> LLH Scans </h2>

Right, we've now pruned our trees and used them to build PDFs and a simulated dataset(s). In theory, we're now ready to launch a fit! But let's be steady Eddies, Cautious Carols, and Nervous Nigels, and make sure everything is behaving as intended. Running fits can be computationally and storage expensive, so it would be a shame to find out after running them that something had previously gone wrong.  

One of way doing this is with likelihood scans. The `llh_scan` app first builds the same test statistic that the fit will use, and we will hand it the same dataset we intend to hand the fitter. For a single floated parameter, it will scan over a range of values centred at the nominal value. For each of 150 steps, it will rebuild the dataset by scaling the PDFs and applying any systematics all set at their nominal values, except the parameter in question which is set to the value we've reached in the scan. This is compared to the actual dataset usinng the test statistic (probably the binned LLH). The LLH is saved at each step, and once we reach the end of the scan for a parameter, it is set back to its nominal value, and we repeat for the next parameter.  

For the 'true' Asimov dataset, these scans should all minimise at 1 (where the x axis is parameter value / nominal value) i.e., by changing a parameter, you can't produce a dataset more similar to the target dataset than by using the exact values used to produce the target dataset. For the other forms of Asimov dataset, most parameters will minimise close to 1, but probably not exactly. The event types with higher rates will be closer to 1 as changes in these will have a greater impact on the LLH, and the small fluctuations causing the differences between the PDFs and Asimov dataset will be smaller proportionally.  

To run the LLH scan:  

>> ./bin/llh_scan ../results/fit_config_file.cfg ../pdfs/config_file.cfg ../cuts/config_file.cfg ../systs/config_file.cfg data_file asimov_rates_file dimensions output_directory

`data_file` is the dataset you want to compare the scaled PDFs to during the scan. This will be the dataset you hand to the fit later. `asimov_rates_file` is a file that contains a mapping of parameter names to Asimov rates i.e., your Asimov dataset file. This is what the PDFs will be scaled by for the nominal values. It will be common for these two files to be the same, however they could be different if for example, you wanted to see how the LLH varied around individual parameters for data. Obviously for data we don't have individual event type rates, so you would give it an Asimov file and scan around those central values.  

`dimensions` is the dimensions of the dataset you want to use when calculating the LLH. This is actually a string and should be "2d" (energy and radius), "3d" (energy, radius, and timePSD), or "4d" (energy, radius, timePSD, and shapePSD). `output_directory` is where you want to save the scan results. A root file `llh_scan.root` will be saved in that directory. In the file will be a plot of LLH vs parameter value for each parameter.  

Running this should be very quick (< 1 minute), though will be longer when using systematics.  

<h2> Running a Fit </h2>

Now the time has come to run a fit. First you want to make sure everything in you fit config is looking sensible. The min and max values, and sigmas for each event type are fairly well tuned, so it's advised you leave these alone unless you know what you're doing. Certainly, for your first fit testing out the code I would leave it as is. If you have good reason to change them, you probably don't need to be reading this guide! `iterations` is the number of steps in the Markov Chain you want to run. `burn_in` is the number of steps you want to cut out the analysis as they occur before the chain has reached the stationary distribution. `fit_distributions` is the event types you want to float in the fit, this will likely be 'all'.  

`output_directory` is where you want to save the fit results, this can be overwritten as an argument to the app. `sigma_scale` is a global value that all parameter's step sizes get scaled by. Again, you probably don't need to change that right now. `hmc_iterations` and `hmc_burn_in` are like `iterations` and `burn_in`, but for running Hamiltonian Markov Chain Monte Carlo. It was found that when running Gaussian step proposals the chain converges more efficiently (more steps but less time per step) than HMCMC (which is generally more efficient for higher numbers of parameters), particularly when including systematics. The functionality to use either or both MCMC and HMCMC was preserved in `fit_dataset`, so you could run some MCMC over a lot of steps quickly, before dialling in with fewer HMCMC steps. This does not seem to be any better than just running MCMC, but the ability to do so remains. Finally, `beeston_barlow` is a bool, when `true` the Barlow Beeston method to account for the MC statistical uncertainty is included in the calculation of the likelihood.  

Once that’s all in line, it may be worth running a small fit with a small number of steps, just to check everything is pointing to the right places. So, set `iterations` and `burn_in` to low numbers, along with their Hamiltonian counterparts. Then run with:  

>> ./bin/fit_dataset ../results/fit_config_file.cfg ../pdfs/config_file.cfg ../cuts/config_file.cfg ../systs/config_file.cfg data_file dimensions output_directory

The config files and `dimensions` are as described for other apps. `output_directory` will override the `output_directory` in the fit config.  

For running full fits, it’s advisable to run multiple fits at once in parallel as you probably want around 1 million steps. Every batch system will be different, but for submitting to the Oxford clusters there is a script `util/submitCondor.py` based on one for submitting RAT jobs originally from Josie (I think).  

You can submit 100 chains by doing something like:  

>> python utils/submitCondor.py dataset_directory dataset_name.h5 fit_name -d 2d -f fit_config_3y.ini -c cut_config.ini -p pdf_config.ini -s syst_config.ini -e /home/parkerw/Software/env-dev.sh -n 100

This will make a directory called `fit_name` within `dataset_directory`. The idea being any different fits you run over the same dataset are all under that directory. Within the `fit_name` directory, there will be subdirectories called `fit_name_i` where `i` runs from 0 to the number of parallel chains you've run (100 in the above example).  

Note that the config filenames here are names and not paths. The python script knows to look in each subdirectory. The environment file you supply here should be whatever you use to set up ROOT, python, GSL, etc., not `bb-likelihood-analysis/env.sh` (which should be sourced before running the python submission.  


<h2> Postfit Analysis </h2>

Now you’re fit has run, the real fun starts! For each individual chain, you’ll have a number of files and subdirectories. In `1dlhproj` and `2dlhproj`, you have projections of the LLH for each parameter, and each combination of two parameters. In each case all other parameters are marginalised over. The burn-in steps are automatically not included.  

`config_log.txt` and `config_log_hmc.txt` save the configs used (these are almost certainly the same as each other). These are also saved in the `cfg` directory if running the python submission script in `utils`. `auto_correlations.txt` calculates how correlated the LLH is to the LLH at a step a different number of steps previous. This should hopefully be close to 0 after a few hundred. `fit_results.txt` contains each parameter value for the maximum LLH step. `scaled_dists` contains the PDFs for each event type, scaled by the event type parameter value at the maximum LLH step in this chain. Also saved is the sum of these, with each systematic’s value at the maximum LLH step applied.  

`error`, `output` and `log` contain the console outputs and batch log. `sh` and `submit` contain the files used to submit these jobs. These all won’t be saved if you ran straight from the terminal (as opposed to using the python submission script in `utils`).  

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
