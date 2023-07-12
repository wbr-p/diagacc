# diagacc

# Overview
Early version of an R Package to estimate diagnostic accuracy.

estimation() is the main function and outputs diagnostic accuracy estimates for the 
area under the curve of the receiver operating characteristic, sensitivity and specificity, 
including confidence intervals.

# Notes
This is an early version of the Package, uploaded in order to be submitted in conjunction with an upcoming publication.

Functionality and convenience features that are planned but not yet implemented include:

1) predictive value estimates
2) compatibility with factorial designs
3) covariate adjustment of estimates
4) compatibility with summary(), plot(), etc.
5) extensive output
6) full documentation

# Installation
The easiest way to install this package is to use the 'devtools' Package from CRAN and execute the command
'install_github("wbr-p/diagacc")'.

# Issues
Should using this Package result in clearly false estimates or you discover bugs of any sort, please file an issue with a minimal reproducible counterexample.
