# ApproximateEM
This repository contins all necessary R scripts to implement the Approximate EM algorithm on joint models of survival and multivariate joint data.
Specifically, this repository fits such models for three, five and ten coviarates each with single time-to-event.

This README briefly sets out file structure. The working directory needs to be set to `~/.../ApproximateEM` in order for the structure set out here to work without a hitch.

Necessary packages are
* `dplyr`
* `lme4`
* `Rcpp` and `RcppArmadillo`
* `survival`
* `ucminf`
* `statmod`

The file `EM.R` fits a trivariate model via approximate EM algorithm. Corresponding programs exist for the five-variate and ten-variate joint models, and are located in the `Kvariate/` folder.

# Data structure
Each longitudinal sub-model is fit with the following formulation:

<img src="https://latex.codecogs.com/svg.latex?{\color{Red}&space;\boldsymbol{Y}_{ik}=(\beta_{k0}&plus;b_{ik0})&plus;(\beta_{k1}&plus;b_{ik1})\times\texttt{time}_{ik}&space;&plus;&space;\beta_{k2}\times\texttt{cont}_i&space;&plus;&space;\beta_{k3}\times\texttt{bin}_i}" title="{\color{Red} \boldsymbol{Y}_{ik}=(\beta_{k0}+b_{ik0})+(\beta_{k1}+b_{ik1})\times\texttt{time}_{ik} + \beta_{k2}\times\texttt{cont}_i + \beta_{k3}\times\texttt{bin}_i}" />

To that end, the algorithm requires data with a column for subject ID called 'id'; time variable 'time'; continuous variable 'cont' and binary 'bin'.  The K longitudinal outcomes should be named 'Y.1, ..., Y.K'. Currently implementation is predicated upon a balanced design only.

# Simulations
Simulations are carried out using the function `simData` `R` package `joineRML` (Hickey *et al.*, 2018). This function generates data in the format as explained above.

The `TrivariateSimulations.R` file generates simulated data using this function under three covariates. These are then fit using the `EM.R` file in the parent folder. The coresponding simulation programs for five- (`K5-`) and ten-covariate (`K10-sims.R`) simulate data for the number of covariates using true values from univariate `lme4` models fit to the ADNI data (see Appendix for further detail). These five and ten-variate version are then fit with corresponding versions of approximate EM found in the `Kvariate/` subdirectory.
