#' Winner's Curse adjustment method using discovery and replication summary statistics - Bowden and Dudbridge (2009)
#'
#' @param summary_disc Data frame containing summary data from discovery GWAS, three columns: rsid, beta, se
#' @param summary_rep Data frame containing summary data from replication GWAS, three columns: rsid, beta, se
#' @param alpha significance threshold used in discovery GWAS
#'
#' @return Data frame with only SNPs deemed significance in discovery GWAS, contains their beta and se values from both discovery and replication GWAS
#' as well as the UMVCUE for each SNP
#' @export
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
