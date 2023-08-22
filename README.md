Simulation Study Comparing Methods for Transportation and Data-Fusion 
=====================================================================

This repository contains code for reproducing the simulation study reported in 'A Calibration Approach to Transportability and Data-Fusion with Observational Data' (Josey et al., 2022). In this numerical experiment, we test and compare different methods for balancing covariate moments between samples and treatments concomitantly estimate the population average treatment effect when there is confounding present. We compare our proposed calibration approaches with a targeted minimum loss approach (Rudolph & van der Laan, 2017) and ab augmented approach that also utilizes calibration weights (Lee et al., 2023).

## R Scripts

- [`main.R`](https://github.com/kevjosey/fusion-sim/tree/master/main.R): Main script for executing simulations. The primary methods we test are all doubly-robust. They include the proposed calibration approach (Josey et al., 2022), an augmented estimator (Lee et al., 2020), and a targeted minimum loss approach (Rudolph & van Der Laan. 2017).
- [`augment.R`](https://github.com/kevjosey/fusion-sim/tree/master/geex.R): Code for fitting doubly-robust estimators described in Lee et al. (2023). Includes a data-fusion and transportability estimator.
- [`tmle.R`](https://github.com/kevjosey/fusion-sim/tree/master/tmle.R): TMLE estimator for the target population average treatment effect from Rudolph & van der Laan (2017).
- [`calibrate.R`](https://github.com/kevjosey/fusion-sim/tree/master/calibrate.R): Methods for fitting balancing weights both in the data-fusion and the transportability cases.
- [`simfun.R`](https://github.com/kevjosey/fusion-sim/tree/master/simfun.R): Contains functions for generating data and for fitting the different estimates. Target and trial samples are generated supposing a sampling score to align with the proposed methodology in the manuscript (Josey et al., 2022). The treatments conditional on the covariates and the sample indicator are constructed using sample specific propensity scores to test the propensity score exchangeability assumption.
- [`unit.R`](https://github.com/kevjosey/fusion-sim/tree/master/unit.R): Unit test the accuracy and functionality of the simulation code.

## References

Lee, D., Yang, S., Dong, L., Wang, X., Zeng, D., & Cai, J. (2023). Improving trial generalizability using observational studies. Biometrics, 79(2), 1213-1225.

Josey, K. P., Yang, F., Ghosh, D., & Raghavan, S. (2022). A calibration approach to transportability and data‐fusion with observational data. Statistics in Medicine, 41(23), 4511-4531.

Rudolph, K.E. and van der Laan, M.J. (2017). Robust estimation of encouragement design intervention effects transported across sites. Journal of the Royal Statistical Society Series B: Statistical Methodology, 79(5), 1509-1525.
