# Winner's Curse Adjustment Methods for GWAS summary statistics


<!-- badges: start -->
  [![R-CMD-check](https://github.com/amandaforde/winners_curse/workflows/R-CMD-check/badge.svg)](https://github.com/amandaforde/winners_curse/actions)
  [![codecov](https://codecov.io/gh/amandaforde/winnerscurse/branch/main/graph/badge.svg?token=P5RVSVWEIM)](https://codecov.io/gh/amandaforde/winnerscurse)
  <!-- badges: end -->



This package has been designed to provide easy access to published methods which aim to correct for Winner's Curse, using GWAS summary statistics. With merely estimates of the association, `beta`, and corresponding standard error, `se`, for each SNP, this package permits users to implement adjustment methods to obtain less biased estimates of the true `beta` values. Methods can be applied to data relating to both quantitative and binary traits.


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

The goal of the functions in this package is to adjust the raw effect estimates, `beta`, rendering them less biased. These adjustments are made using only the summary statistics obtained from the discovery study.  



### Simple Example

We will now discuss a simple simulated data set, providing another way of comprehending the Winner's Curse with the assistance of a plot. However, remember throughout that this is very much a hypothetical data set used for demonstration purposes and does not resemble a true data set that one would obtain in a GWAS.  

Let us first generate 1000 values from a standard normal distribution, centred at 0 with standard deviation of 1. We allow these values to act as the *true* effect sizes of 1000 independent SNPs.

```r
set.seed(1948)
beta <- stats::rnorm(1000,0,1)
beta <- data.frame(id = 1:1000, beta = beta)
````

Next, we use these true effect sizes to simulate *one* effect estimate for each SNP - this is similar to performing a single GWAS study or looking at points scored by players in one World Cup tournament. We imagine that the effect size of each SNP also follows a normal distribution which is centred at its true effect size, `beta` and has a standard deviation of 0.5. This value of 0.5 is only used for convenience and in general, each SNP would have its own individual value here. Therefore, `beta_hat` below can be understood as a vector of estimates of effect sizes obtained from one study. We then re-arrange our data set so that SNPs are now ordered based on decreasing absolute `beta_hat` values and have a look at the 6 most extreme SNPs. As this data set is simulated, we know the true `beta` values and thus, can easily compare `beta_hat` with `beta` for each SNP. With real data sets, we are not as fortunate. From the table below, the bias induced by Winner's Curse is evident. Each `beta_hat` value is an *overestimation* of its corresponding true effect size.  

``` r
beta_hat <- stats::rnorm(1000,beta$beta,0.5)  
data <- cbind(beta, beta_hat)
data <- dplyr::arrange(data, desc(abs(beta_hat)))
head(data)
```


![](https://raw.githubusercontent.com/amandaforde/winnerscurse/main/winnerscurse.PNG)



Plotting the distribution of the true absolute effect sizes for the SNPs which have been ranked in the top 50 together with the estimated absolute effect sizes from this study allows us visualize Winner's Curse. As expected, we witness the orange curve of the estimated effects lying to the right of the true effect turquoise curve. This provides a strong indication that a large number of the estimated effects of these extreme SNPs are greater than their corresponding true value. 

```r
library(RColorBrewer)
col <- brewer.pal(3,"Dark2")
plot(density(abs(data$beta[1:50])),ylim=c(0,1.4),xlim=c(0.5,4.5), main="Visualisation of Winner's Curse", col=col[1], lwd=2)
lines(density(abs(data$beta_hat[1:50])), col=col[2], lwd=2)
```


![](https://raw.githubusercontent.com/amandaforde/winnerscurse/main/readme_plot.png)


Therefore, in visual terms, we would like the functions of our package to produce estimates which are more in line with the true values, i.e. shifted more towards the turquoise density plot on the left, avoiding the obvious inflation incurred by the raw `beta_hat` estimates here. 

**Note:** The above data is *not* suitable for use with the functions in this package. It has merely been used a very simple example to demonstrate the concept of winner's curse. In order to appropriately use the functions, summary statistics must be in the form of a **data frame** in which the first column, titled `rsid`, contains the SNP ID number, the second column, named `beta`, contains the effect size estimate while the third column, `se` holds the corresponding estimated standard error. 
