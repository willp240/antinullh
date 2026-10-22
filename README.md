Reactor Antinu OXO Analysis
=======================

This repository contains code to run a likelihood analysis of reactor antineutrinos at SNO+. The fitting is done using `OXO`.

The four main executables use two models:

| Executable | Model | Purpose |
| --- | --- | --- |
| `mcmc` | Floating oscillation through a shape systematic | Sample the posterior with Metropolis MCMC |
| `llh_scan` | Same model as `mcmc` | Scan one parameter at a time with the others nominal |
| `fixedosc_fit` | Oscillation fixed in the reactor PDFs | Fit the remaining parameters with Minuit |
| `fixedosc_llhscan` | Same model as `fixedosc_fit` | Scan one parameter at a time, rebuilding reactor PDFs for oscillation points |

Both scans show likelihood slices through the nominal point, not profile likelihoods. Shared analysis code is under `src/`; each executable's run-specific logic is under `exec/`.

<h3>Installation:</h3>

Clone this repository \
Edit the env.sh file to point to your `OXO` install, this `antinullh` directory, and the top level directory where you want to files to be written to

> source env.sh \
> make

<h3>Running:</h3>

Full instructions on how to run the code is in the `HOWTO.md` guide
