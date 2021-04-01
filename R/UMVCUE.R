#'UMVCUE method for use with discovery and replication GWASs
#'
#'\code{UMVCUE} is a function which aims to produce less biased SNP-trait
#'association estimates for SNPs deemed significant in the discovery GWAS, using
#'summary statistics from both discovery and replication GWASs. The function
#'implements the method described in
#'\href{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2726957/}{Bowden and
#'Dudbridge (2009)}, which was established for this purpose.
#'
#'
#'@param summary_disc A data frame containing summary statistics from the
#'  \emph{discovery} GWAS. It must have three columns with column names \code{rsid},
#'  \code{beta} and \code{se}, respectively, and all columns must contain
#'  numerical values.
#'@param summary_rep A data frame containing summary statistics from the
#'  \emph{replication} GWAS. It must have three columns with column names \code{rsid},
#'  \code{beta} and \code{se}, respectively, and all columns must contain
#'  numerical values.
#'@param alpha A numerical value which specifies the desired genome-wide
#'  significance threshold for the discovery GWAS.
#'
#'@return A data frame with summary statistics and adjusted association
#'  estimate of only those SNPs which have been deemed significant in the
#'  discovery GWAS according to the specified threshold, \code{alpha}, i.e. SNPs
#'  with \eqn{p}-values less than \code{alpha}. The inputted summary data
#'  occupies the first five columns, in which the columns \code{beta_disc} and
#'  \code{se_disc} contain the statistics from the discovery GWAS and columns
#'  \code{beta_rep} and \code{se_rep} hold the replication GWAS statistics. The
#'  new adjusted association estimate for each SNP, as defined in the
#'  aforementioned paper, is contained in the final column, namely
#'  \code{beta_UMVCUE}. The SNPs are contained in this data frame according to
#'  their significance, with the most significant SNP, i.e. the SNP with the
#'  largest absolute \eqn{z}-statistic, now located in the first row of the data
#'  frame.
#'@references Bowden, J., & Dudbridge, F. (2009). Unbiased estimation of odds
#'  ratios: combining genomewide association scans with replication studies.
#'  \emph{Genetic epidemiology}, \strong{33(5)}, 406\eqn{-}418.
#'  \url{https://doi.org/10.1002/gepi.20394}
#'
#'@seealso
#'
#'\url{https://amandaforde.github.io/winnerscurse/articles/discovery_replication.html}
#'for illustration of the use of \code{UMVCUE} with toy data sets and further
#'information regarding computation of the adjusted SNP-trait association
#'estimates for significant SNPs.
#'@export
#'
#'
UMVCUE <- function(summary_disc, summary_rep, alpha){

  c <- stats::qnorm(1-alpha/2)

  disc_sig <- summary_disc[abs(summary_disc$beta/summary_disc$se) > c,]
  rep_sig <- summary_rep[abs(summary_disc$beta/summary_disc$se) > c,]

  disc_z <- abs(disc_sig$beta/disc_sig$se)
  disc_sig <- cbind(disc_sig,disc_z)
  rep_sig <- cbind(rep_sig,disc_z)

  disc_sig <- dplyr::arrange(disc_sig,dplyr::desc(disc_sig$disc_z))
  rep_sig <- dplyr::arrange(rep_sig,dplyr::desc(rep_sig$disc_z))


  mu_hat <- function(i){return((((rep_sig$se[i])^2)*(disc_sig$beta[i]) + ((disc_sig$se[i])^2)*(rep_sig$beta[i]))/((disc_sig$se[i])^2 + (rep_sig$se[i])^2))}
  w <- function(s,t,p){return(((sqrt((disc_sig$se[s])^2 + (rep_sig$se[s])^2))/((disc_sig$se[s])^2))*(mu_hat(s) -((-1)^p)*((disc_sig$se[s]*abs(disc_sig$beta[t]))/(disc_sig$se[t]))))}
  beta_UMVCUE_i <- function(i){return(mu_hat(i) - (((rep_sig$se[i])^2)/(sqrt((disc_sig$se[i])^2 + (rep_sig$se[i])^2)))*((stats::dnorm(w(i,i+1,0))-stats::dnorm(w(i,i-1,0))-stats::dnorm(w(i,i+1,1))+stats::dnorm(w(i,i-1,1)))/(stats::pnorm(w(i,i+1,0))-stats::pnorm(w(i,i-1,0))-stats::pnorm(w(i,i+1,1))+stats::pnorm(w(i,i-1,1)))))}
  w_0 <- function(s,t,p){return(((sqrt((disc_sig$se[s])^2 + (rep_sig$se[s])^2))/((disc_sig$se[s])^2))*(mu_hat(s) -((-1)^p)*((disc_sig$se[s]*abs(Inf)))))}
  beta_UMVCUE_1 <- function(i){return(mu_hat(i) - (((rep_sig$se[i])^2)/(sqrt((disc_sig$se[i])^2 + (rep_sig$se[i])^2)))*((stats::dnorm(w(i,i+1,0)) - stats::dnorm(w_0(i,i-1,0))-stats::dnorm(w(i,i+1,1))+stats::dnorm(w_0(i,i-1,1)))/(stats::pnorm(w(i,i+1,0))-stats::pnorm(w_0(i,i-1,0))-stats::pnorm(w(i,i+1,1))+stats::pnorm(w_0(i,i-1,1)))))}
  beta_UMVCUE_N <- function(i){return(mu_hat(i) - (((rep_sig$se[i])^2)/(sqrt((disc_sig$se[i])^2 + (rep_sig$se[i])^2)))*((stats::dnorm(w_0(i,i+1,0)) - stats::dnorm(w(i,i-1,0))-stats::dnorm(w_0(i,i+1,1))+stats::dnorm(w(i,i-1,1)))/(stats::pnorm(w_0(i,i+1,0))-stats::pnorm(w(i,i-1,0))-stats::pnorm(w_0(i,i+1,1))+stats::pnorm(w(i,i-1,1)))))}

  n_sig <- length(disc_z)
  beta_UMVCUE <- c(rep(0,n_sig))

  beta_UMVCUE[1] <- beta_UMVCUE_1(1)
  beta_UMVCUE[n_sig] <- beta_UMVCUE_N(n_sig)


  for (i in 2:(n_sig-1)){
    beta_UMVCUE[i] <- beta_UMVCUE_i(i)
  }

  summary_data_sig <- cbind(disc_sig[1:3],rep_sig[2:3],beta_UMVCUE)
  names(summary_data_sig)[2] <- "beta_disc"
  names(summary_data_sig)[3] <- "se_disc"
  names(summary_data_sig)[4] <- "beta_rep"
  names(summary_data_sig)[5] <- "se_rep"

  return(summary_data_sig)

}
