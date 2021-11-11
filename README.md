# ApproximateEM
This repository contins all necessary R scripts to implement the Approximate EM algorithm on joint models of survival and multivariate joint data.

This README briefly sets out file structure. The working directory needs to be set to `~/.../ApproximateEM` in order for the structure set out here to work without a hitch.

Necessary packages are
* `dplyr`
* `nlme`
* `Rcpp` and `RcppArmadillo`
* `survival`
* `ucminf`
* `statmod`

# Function call and output

The file `EM.R` fits a multivariate model with K longitudinal responses via approximate EM algorithm. The user can determine the number of responses via the `nK` argument in the call to `em()`.

## User input
Necessary arguments here are:
* `data`: a dataset that matches the required data format detailed in the next section. Currently only completely balanced data is fit by the algorithm;
* `ph`: a Cox proportional-hazards model fit with a continuous `cont` and binary `bin` variable which are additionally included in the longitudinal model and data;
* `nK`: The number of longitudinal responses to fit. (if nK=3 then variables named Y.1 -> Y.3 are fit in a trivariate joint model.
Optional arguments are:
* `gh.nodes`: The number of weights and abscissae to use in Gauss-Hermite quadrature (defaults to three);
* `collect.hist`: Should parameter estimates at each iteration be stored? (defaults to `TRUE`);
* `max.iter`: How many iterations to attempt before exiting and returning the max.iter'th as final parameter estimates (defaults to 200);
* `tol`: Tolerance that difference between iterations across all parameters must be less than (defaults to 0.01);
* `diff.type`: Should relative difference (the defualt, 'abs.rel') be used, or absolute ('abs') be used to monitor convergence?
* `post.process`: Should posterior estimates of random effects be calculated at final parameter estimates and standard errors calculated? (defaults to `TRUE`);
* `verbose`: If `TRUE` then prints all parameter estimates to the console at each iteration (defaults to `FALSE`).

## Generated Output
* `REs`: A matrix of the random effects. If `post.process` is `TRUE` then these are calculated post-hoc at final parameter values;
* `coeffs`: A list of parameter estimates;
* `EMtime`: The time taken for the approximate EM algorithm to converge (s);
* `mvlme.time`: The time taken for the MVLME step to obtain initial conditions to converge (s);
* `comp.time`: The time taken for the `em()` call to complete (i.e. EM, MVLME, data functions, all initial conditions and post-processing;
* `history`: if `collect.hist` is `TRUE` then a matrix containing parameter estimates at each iteration.
* `SEs`: if `post.process` is `TRUE` then a vector of approximate standard errors calculated using the observed empirical information matrix.
* `postprocess.time`: if `post.process` is `TRUE` then the time taken for RE calculation and standard error calculation in post processing step.

# Data structure
Each longitudinal sub-model is fit with the following formulation:

<img src="https://latex.codecogs.com/svg.latex?{\color{Red}&space;\boldsymbol{Y}_{ik}=(\beta_{k0}&plus;b_{ik0})&plus;(\beta_{k1}&plus;b_{ik1})\times\texttt{time}_{ik}&space;&plus;&space;\beta_{k2}\times\texttt{cont}_i&space;&plus;&space;\beta_{k3}\times\texttt{bin}_i}" title="{\color{Red} \boldsymbol{Y}_{ik}=(\beta_{k0}+b_{ik0})+(\beta_{k1}+b_{ik1})\times\texttt{time}_{ik} + \beta_{k2}\times\texttt{cont}_i + \beta_{k3}\times\texttt{bin}_i}" />

To that end, the algorithm requires data with a column for subject ID called 'id'; time variable 'time'; continuous variable 'cont' and binary 'bin'. The K longitudinal outcomes should be named 'Y.1, ..., Y.K'. Currently implementation is predicated upon a balanced design only.

# Simulations
Simulations are carried out using the function `simData` `R` package `joineRML` (Hickey *et al.*, 2018). This function generates data in the format as explained above.

The `TrivariateSimulations.R` file generates simulated data using this function under three covariates. These are then fit using the `EM.R` file in the parent folder. The coresponding simulation pro for five (`K5-`) and ten-covariate (`K10-sims.R`) simulate data for the number of covariates using true values from univariate `lme4` models fit to the ADNI data (see Supplementary Materials for further detail). These five and ten-variate version are then fit using the `em()` function with argument `nK` set to 5 or 10.
