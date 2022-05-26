<!-- badges: start -->
  [![R-CMD-check](https://github.com/amandaforde/winners_curse/workflows/R-CMD-check/badge.svg)](https://github.com/amandaforde/winners_curse/actions)
  [![codecov](https://codecov.io/gh/amandaforde/winnerscurse/branch/main/graph/badge.svg?token=BORRC1EUZ7)](https://codecov.io/gh/amandaforde/winnerscurse)
<!-- badges: end -->


# Winner's Curse Adjustment Methods for GWAS summary statistics

This package, `winnerscurse`, has been designed to provide easy access to published methods which aim to correct for Winner's Curse, using GWAS summary statistics. With merely estimates of the SNP-trait association, `beta`, and corresponding standard error, `se`, for each SNP, this package permits users to implement adjustment methods to obtain less biased estimates of the true `beta` values. Methods can be applied to data relating to both quantitative and binary traits. This package contains functions which can be implemented with just the summary statistics from a discovery GWAS as well as functions which require summary data from both discovery and replication GWASs. Users can also obtain confidence intervals and standard errors for certain adjusted association estimates. 


### Installation

You can install the current version of `winnerscurse` from GitHub with:

```r
# install.packages("devtools")
devtools::install_github("amandaforde/winnerscurse")
```


### Winner's Curse 

The Winner's Curse is a statistical effect usually resulting in the exaggeration of SNP-trait association estimates in the sample in which these associations were discovered. 

However, the Winner's Curse isn't just a phenomenon related to GWAS. An understanding of Winner's Curse may be gained through the following simple example. Consider rugby players who are ranked based on the number of points scored in one World Cup tournament. The top ranking players, the *'winners'*, are most likely players who had an above average tournament, perhaps those who were at the peak of their career and fitness. However, across many tournaments, these highly ranked players may not score as many points nor be as outstanding consistently. 

Now, let us reframe this idea in the context of GWAS. The *'winners'* here are SNPs whose effect sizes are *stochastically* higher in the discovery study than their true association values. Clearly, these raw effect estimates are therefore biased estimates of `beta`, especially for those SNPs who are ranked highly in the study - the most significant SNPs. SNPs are often ranked according to their *z*-statistics or corresponding *p*-values. 

The goal of the functions in this package is to adjust the raw effect estimates, `beta`, rendering them less biased. Several functions in this package can make adjustments using only the summary statistics obtained from the discovery study.  



### Example

We will now discuss a simple simulated data set of GWAS summary statistics, providing an alternative manner in which to comprehend the Winner’s Curse with a visual representation of the bias induced. We simulate GWAS summary statistics for a quantitative trait with a normal effect size distribution under the assumption that SNPs are independent. We assume a fixed array of 1,000,000 SNPs, and fix the heritability of this trait at 0.7, the sample size to be 30,000 and the proportion of effect SNPs at 0.01. A more detailed description of the code used to simulate this data set can be viewed in the article [‘Methods for use with discovery GWAS’](https://amandaforde.github.io/winnerscurse/articles/winners_curse_methods.html). 
The true effect sizes, `beta`, were first generated and then, used these to simulate *one* effect estimate, `beta_hat` for each SNP – this is similar to performing a single GWAS or only having access to the points scored by players in one World Cup tournament. In this simulation scenario, we can thus easily compare `beta_hat` with `beta` for each SNP. With real data sets, we are not as fortunate and hence, this provides motivation to establish suitable Winner’s Curse correction methods which will aim to reduce or eliminate the bias witnessed in the effect estimates of significant SNPs. The table below shows the 6 SNPs which have been deemed most extreme according to their z-statistic. It contains the values for `beta_hat` and `beta`, with the bias induced by Winner’s Curse evident. 5 of the `beta_hat` values are an *overestimation* of their corresponding true effect size.

<p align="center">
  <img src="https://raw.githubusercontent.com/amandaforde/winnerscurse/main/winnerscurse_table.PNG" width="25%">
</p>

Plotting the distribution of the true absolute effect sizes for SNPs which have been deemed significant according to the genome-wide significance threshold of `5e-8` with the estimated absolute effect sizes from this study allows visualise Winner’s Curse. As expected, we witness the orange curve of the estimated effect lying to the right of the true effect turquoise curve. This provides a strong indication that a large number of the estimated effect of these extreme SNPs are greater than their corresponding true value. 

<p align="center">
  <img src="https://raw.githubusercontent.com/amandaforde/winnerscurse/main/winnerscurse_plot.PNG" width="50%">
</p>

Therefore, in visual terms, we would like the functions of our package to produce estimates which are more in line with the true values, i.e. shifted more towards the turquoise density plot on the left, avoiding the obvious inflation incurred by the raw `beta_hat` estimates here. 

**Note:** In order to appropriately use the functions in this package, summary statistics must be in the form of a **data frame** in which the first column, titled `rsid`, contains the SNP ID number, the second column, named `beta`, contains the effect size estimate while the third column, `se` holds the corresponding estimated standard error. 
