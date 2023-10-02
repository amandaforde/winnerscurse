#' Simulating GWAS summary statistics
#'
#' \code{sim_stats} is a function which can be used to simulate summary
#' statistics for a set of independent SNPs for both discovery and replication
#' GWASs. This function allows the user to create toy datasets with which they
#' can explore the implementation of the Winner's Curse correction methods.
#'
#' @param nsnp A numerical value which specifies the total number of SNPs that
#'   the user wishes to simulate summary statistics for. The default is 1,000,000 SNPs,
#'   i.e. \code{nsnp=10^6}.
#' @param h2 A numerical value between 0 and 1 which represents the desired
#'   heritability of the trait of interest, or in other words, the total
#'   variance explained in the trait by all SNPs. The default is a moderate
#'   heritability value of 0.4, \code{h2=0.4}.
#' @param prop_effect A numerical value between 0 and 1 which determines the
#'   trait's polygenicity, the fraction of the total number of SNPs which are truly associated with the
#'   trait. The default setting is \code{prop_effect = 0.01}.
#' @param nid A numerical value which specifies the number of individuals that
#'   the discovery GWAS has been performed with. This value defaults to 50,000 individuals, \code{nid=50000}.
#' @param rep A logical value which allows the user to state whether they would
#'   also like to simulate summary statistics for a replication GWAS based on
#'   the same parameters and true effect sizes. The default setting is
#'   \code{rep=FALSE}.
#' @param rep_nid A numerical value which specifies the number of individuals
#'   that the replication GWAS has been performed with. Similar to \code{nid},
#'   this value defaults to 50,000 individuals, \code{nid=50000}.
#'
#' @return A list containing three different components in the form of data
#'   frames, \code{true}, \code{disc} and \code{rep}. The first element,
#'   \code{true} has two columns, \code{rsid} which contains identification numbers
#'   for each SNP and \code{true_beta} which is each SNP's simulated true
#'   association value. \code{disc} has three columns representing the summary statistics
#'   one would obtain in a discovery GWAS. For each SNP, this data frame contains
#'   its ID number, its estimated effect size, \code{beta}, and associated standard error, \code{se}.
#'  Similarly, if the \code{rep} argument in the function has been set to \code{TRUE},
#'  then the data frame, \code{rep} has three columns representing the summary statistics
#'   one would obtain in a replication GWAS. In this data frame, just as with \code{disc},
#'   the values for \code{beta} have been simulated using the true association values, \code{true_beta},
#'   and the standard errors are reflective of the chosen sample size.
#'   If \code{rep=FALSE}, \code{NULL} is merely returned for this third element.
#' @export
#'


sim_stats <- function(nsnp=10^6,h2=0.4,prop_effect=0.01,nid=50000,rep=FALSE,rep_nid=50000){
  effect_snps <- prop_effect*nsnp
  maf <- stats::runif(nsnp,0.01,0.5)
  true_beta <- stats::rnorm(effect_snps,0,sd=sqrt((2*maf*(1-maf))))
  var_y <- sum(2*maf[1:effect_snps]*(1-maf[1:effect_snps])*true_beta^2)/h2
  true_beta <- true_beta/sqrt(var_y)
  true_beta <- c(true_beta, rep(0,nsnp-effect_snps))
  se <- sqrt((1 - 2*maf*(1-maf)*true_beta^2)/(2*(nid-2)*maf*(1-maf)))

  true <- data.frame(rsid=seq(1,nsnp),true_beta=true_beta)

  disc <- data.frame(rsid=seq(1,nsnp),beta=stats::rnorm(n=nsnp,mean=true_beta,sd=se),se=se)

  if(rep==TRUE){
    se_rep <-  sqrt((1 - 2*maf*(1-maf)*true_beta^2)/(2*(rep_nid-2)*maf*(1-maf)))
    rep <- data.frame(rsid=seq(1,nsnp),beta=stats::rnorm(n=nsnp,mean=true_beta,sd=se_rep),se=se_rep)
  }else{rep <- NULL}

  return(list("true"=true,"disc"=disc,"rep"=rep))
}
