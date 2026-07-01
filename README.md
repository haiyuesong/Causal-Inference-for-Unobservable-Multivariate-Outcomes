# Causal Inference for Unobservable Multivariate Outcomes, with Applications to Brain Effective Connectivity
## Information
- Authors: Haiyue Song, Ani Eloyan, Youjin Lee


## Code
**Simulation studies** for evaluating finite-sample performance under null and alternative settings. 
- `simulation studies/ec_kernals.cpp`: C++ implementation of core computational functions used in the simulation studies, including efficient calculations related to effective connectivity estimation.
- `simulation studies/simulation_null.R`: Runs simulation under the global null hypothesis, primarily used to assess type-I error and familywise error control.
- `simulation studies/simulation_alternative.R`: Runs simulation under alternative hypotheses, used to evaluate false discovery proportion and power in detecting nonzero causal effects.
- `simulation studies/visualization.R`: Generates figures (Figures 3 and A1) and table (Table A1) for the simulation results.

**Data application** code for preprocessing resting-state fMRI data, computing effective connectivity, constructing the analysis cohort, and conducting causal analyses. 
- `data application/fMRI_preprocessing.R`: Preprocesses resting-state fMRI data prior to effective connectivity analysis.
- `data application/effective_connectivity_computation.R`: Computes subject-level brain effective connectivity measures from preprocessed fMRI time series.
- `data application/functions.R`: Contains helper functions used throughout the data application pipeline.
- `data application/cohort.R`: Constructs the analytic cohort and prepares subject-level covariates and study variables.
- `data application/causal_analysis.R`: Performs the causal inference analysis using the derived effective connectivity outcomes and visualizes the results (Figure 4). Implements diagnostics and sensitivity analysis. 

## Citation
Song, H., Eloyan, A., & Lee, Y. (2026). Causal Inference for Unobservable Multivariate Outcomes, with Applications to Brain Effective Connectivity. arXiv:2604.00390.

```bibtex
@article{song2026causal,
  title   = {Causal Inference for Unobservable Multivariate Outcomes, with Applications to Brain Effective Connectivity},
  author  = {Song, Haiyue and Eloyan, Ani and Lee, Youjin},
  year    = {2026},
  eprint  = {2604.00390},
  archivePrefix = {arXiv},
  primaryClass  = {stat.ME},
  url     = {https://arxiv.org/abs/2604.00390}
}
```
