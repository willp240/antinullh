<h1>Antinu Analysis Event Recon How To Guide</h1>
<h2>Introduction</h2>
This is a preliminary version of SNO+ 2p2 antinu reconstruction analysis. In general, the goal of this analysis is to construct prompt energy corrected spectrum before being used in following oxo fit to get systematics and statistical uncertainty. Analysis is constituted by four parts: 

- tagging 
- Accidental/IBD Classifier
- Efficiency Study and Expected Events Calculation
- Sideband Study

<h2>Version of mc/data</h2>

The data and mc usage is from 2p2 fulfill antinulist. 
In simulation, relevant backgrounds are `Geoibd-U`, `Geoibd-Th` and `Alphan_Lab_13c`, `accidentals` and `BiPo214`

<h2>Taggings</h2>

Main codes are in `powei/Tagging`. Firstly, both data and mc have to pass selection mainly developed from Anthony Zummo. The command is as following:

> path/to/antinurecon.exe path/to/rawinput.root runid path/to/outputfilename.root ratversion 0/1

In the last bit, 0 is for sim whereas 1 is for data. Outputs of this operation also depends on assignation of sim or data. In sim, only `.root` will be produced but a `.root` and a `.txt` recording infor of all >3000 events will be produced if user runs in data mode. 
Another tagging code, `IBDPrompt_DelayRecon.exe`, is for producing accidental source. Note that input root file of this execution has to be assigned as a single subrun, instead of parsing a run like `antinurecon.exe`. Usage of this file is
> path/to/IBDPrompt_DelayRecon.exe path/to/rawinput_subrun.root runid path/to/outputfilenam.roote ratversion 1

A reconstructed `.root` file for that specified subrun will be produced after this command.

<h2>Accidental Study</h2>

Coincidence Tagging more or less produces accidental events. In antinu, this event is not a major background contribution but a classifier can be built from discrepancy of `delayedEcorr` and `dtdr` distribution of reactoribd and accidentals to reduce accidentals being a negligible backgrounds in final spectrum.

<h3>Accidentals Generation</h3>

The first step is to generate millions of accidentals events using data-driven method. In other words, previous prompt and delay events from IBDPrompt_DelayRecon.exe act as a input source of the generator. To run a accidental generator:

> path/to/AccidentalGenerator.cpp path/to/gridprocessantinulist.dat runid path/to/outputfilename.root ratversion 10000

Two things worth mentioning here: 1. The first argument is a `.dat` from `rat-tools/GridTools/processing_list` acting on 2p2antinulist. We need this to point a correct subrun data and starting randomly sampling. 2. The last argument could be any integer number, This is a number to tell generator approximately how many accidentals you want to generate.
<h3>Produce Accidentals' and Reactoribds'(Unoscillated) delayedEcorr and drdt pdfs</h3>

Two programs are in charge of building delayedEcorr and drdt separately.
In overall, one of their inputs is those generated events filepath or unoscillated reactoribd path, in order to build correct shape of distributions. To tun these codes:

> python3 Acc_Reac_delayedE.py
> root -l finalAcc_Reac_dtdR.cpp

**Note that user has to change unoscillated reactoribd and accidentals filepath to their own directories.** A set of `.png/.pdf` and a `.root` with TH1D/TH2D objects are outputted from these programs.

<h3>Calculate IBD/Acc likelihood ratio and Attach new branches into original ntuple files</h3>

Based on distributions above, a likelihood ratio between accidentals and reactoribd could be calculated. This likelihood could be utilised to all interesting eventtypes in simulation, data and accidentals to perform efficiencies studies on these eventtypes later.  

> path/to/ApplyPosterior_Distr.exe path/to/inputntuple.root(aftercuts) path/to/outputntuple.root ratver <data/Acc/Evttype>

**Note that user has to change filepath to their own directories.**
User could pick data, Acc or interesting Evttype such as Geoibd-Th as the last part of arguments. The program will get filepath from given eventtype and save `Log_Llh_IBD` `Log_Llh_Acc` and `post_prob` branches into ntuples.

<h3>Calculate FOM</h3>

A FOM is performed to get the best cut of `post_prob`.

> python3 `ClassifierFOM.py`

It will make a plot of FOM distribution based on a grid-search. Also, It will print out the best cut from FOM,and showing Posterior probability( aka classifier value) for both Reactoribd and Accidentals as a check. 




<h2>Efficiency Study and Expected Events Calculation</h2>

<h3>Tagging Efficiency Study</h3>
To estimate how many events of each eenttypes into spectrum, it is important to do a efficiency study. `effcalculation.cpp` is used to copt with this. User can alter the cuts inside to check efficiencies on cuts applied. However, officially, the tagging efficiency is defined as num of event tagged in 5.7FV/num of simulated = num of event tagged in 5.7FV/eVIndex == 0. 
To check efficiency of a specifice evttype:

> root -l 'effcalculation.cpp(708,"evttype")'
**Note that user has to change filepath to their own directories of rawntuple and reconstructed one.**
Therefore, every evettypes' efficiency should be investigated completely.


Logically, Expected amounts of events are calculated as num of events reconstructed in sim/30000. * livetime in that run, and summing over all runs in antinulist. Therefore, getting livetime per run is crucial to estimate expected event numbers. 
<h3>Livetime Calculation for each run</h3>

We attempt to use official livetime calculator written from Max Smiley.  
<h4>Generate files with a single runid in each from antinulist (Accidentally Deleted ><)</h4>

Firstly, a python script aims to produce thousands of `.txt` files. The only content in each file is a runid. These file will be ready to be parsed as one of the arguments in `scintillator_livetime_calculator.py`.

<h4>Producing livetime for each run</h4>
A command is introduced to get livetime per rum:

> python path/to/scintillator_livetime_calculator.py -i path/to/singlerunidfile -o outputdirpath -d /path/to/rawinputfilepath --veto_dir path/to/hinhitsveto.txt

This program will swallow a singlerunidfile to decide which run it is ready to look at. Then a 20s muon veto will be automatically done inside this program and combine livetime lost with a hinhitsveto file we acquired in tagging (look again in Taggings section if get confused). With all infomation above, several `.txt` tables are outputted but we only care about `livetime_total_defaultlivetime.txt`. 



<h3>Expected Event Numbers Calculation</h3>

Expected events of `Geoibd-Th`, `Geoibd-U` and `Reactoribd` are calculated here. 
**Note that user has to change filepath to their own directories of reconstructed ntuple and livetimepath.**
> python3 ExpectedNumEstimation.py

The program needs to take livetime and amounts of simulation ntuple per run to do countings. Final Expected numbers will be printed out in terminal once user types this command.


<h2> Sideband Study </h2>

To validate our distribution is trustable and demonstrated that we've truly understood every bit of our spectrum, a sideband study is perform in antinu analysis. Sideband region is defined as [0.6, 1.85] MeV, as there is no clear excess above 2.5 MeV.
Before a (\alpha,p) generator is established, sideband is studied as the following way:
Hypothesis is that the excess event between 1.2 and 1.85 might be induced from Po214 alpha interation. To fit this contribution, a tail is "borrowed" from Po215 alpha and shifted linearly by comparing two alphas' Q-value is attached to current Po214 alpha model. Therefore, ingredients of the sideband fit are:

- `Po215 alpha`
- `Po212 alpha`
- `Po214 alpha`
- `Po214 alpha+gamma`
- `Po214 Tail from 215 prefit`
- `Neutron Capture (NC)`
The data is from using antinucuts but freeing low delayedEcorr cut, without BiPo214 veto and Acc/IBD classifier. Here a step-by-step instructions are given to reach the final fit

<h3> Reconstruct Rn219-Po215 </h3>

To gain more statistics in the tail region of Po215, we included partialfill data in antinupartialfilllist to do tagging as well. Two programs are dealing with tagging for 2p2fulfill and partialfill separately:

> /path/to/newRn219Po215Recon.exe path/to/rawinput.root runid path/to/outputfilename.root rat-7.0.8 1
> /path/to/partialRn219Po215Recon.exe path/to/rawinput.root runid path/to/outputfilename.root rat-6.18.9 1

`.root` ntuple file will be yielded.

<h3> Bg Simulation </h3>

Unlike Po215 and NC event using data-driven method, rest of backgrounds are fitted with simulated delayedEcorr distribution. User can do this part from themselves, but one trick is given here for simulating Po214 alpha+gamma as the rate is too low to build Asimov dataset. This credit should give to Valentina. She suggested that user could temporally turn off/delete alph decay line from `beta_decays.dat` in rat. In the meantime, a change of macro is needed to push rat generator follow decaychain methode to do simulation: Turn of decay0 generator in ordinary Po214 generation macro and change to 

```
/generator/add decaychain 214Po:fill:poisson:alpha
/generator/pos/set 0 0 0
``` 
<h3> Fit each Bg distributions respectively</h3>

- `Po215 alpha` 
> root -l sideband_copy.cpp
- `Po215 Tail`
> root -l Sideband_ROI_KDEFit.cpp
- `Po212 alpha`
> root -l Po212alphaFit.cpp
- `Po214 alpha`
> root -l Po214alphaFit.cpp
- `Po214 alpha+gamma`
> root -l Po214alphagammaFit.cpp
- `Neutron Capture (NC)`
- `root -l delayEGausFits.cpp`
Each fits produces a `RooWorkspace` in a `.root` file.


<h3> Final Fit</h3>

Again, the final fit is fitting with a spectrum selected from antinucuts but freeing low delayedEcorr cut, BiPo214 veto and Acc/IBD classifier. It uses those pdfs mentioned above to look at Po214 contribution into ROI. To do this:

> root -l AntiCut214Tag_DelayESpectrumFitUpdate.cpp


