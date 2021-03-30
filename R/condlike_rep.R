#' Conditional likelihood method for Winner's Curse adjustment using discovery and replication summary statistics - Zhong and Prentice (2008)
#'
#' @param summary_disc Data frame containing summary data from discovery GWAS, three columns: rsid, beta, se
#' @param summary_rep Data frame containing summary data from replication GWAS, three columns: rsid, beta, se
#' @param alpha significance threshold used in discovery GWAS
#'
#' @return Data frame with only SNPs deemed significance in discovery GWAS, contains their beta and se values from both discovery and replication GWAS
#' as well as the combined unadjusted estimate, the conditional likelihood MLE and the weighted average for each SNP
#' @export
#'
condlike_rep <- function(summary_disc,summary_rep,alpha){
  ### only selection on snps from discovery data set, just require definition for c_1 or equivalently, alpha

  se_com <- sqrt((((summary_disc$se)^2)*((summary_rep$se)^2))/(((summary_disc$se)^2) + ((summary_rep$se)^2)))
  beta_com <- ((((summary_rep$se)^2)*(summary_disc$beta))+(((summary_disc$se)^2)*(summary_rep$beta)))/(((summary_disc$se)^2) + ((summary_rep$se)^2))
  summary_com <- data.frame(rsid=summary_disc$rsid,beta=beta_com,se=se_com)
  c_1 <- stats::qnorm(1-(alpha)/2)


  summary_disc_sig <- summary_disc[abs(summary_disc$beta/summary_disc$se) > c_1,]
  summary_rep_sig <- summary_rep[abs(summary_disc$beta/summary_disc$se) > c_1,]
  summary_com_sig <- summary_com[abs(summary_disc$beta/summary_disc$se) > c_1,]

  beta_adj <- function(i){
    beta_com <- summary_com_sig$beta[i]
    beta_1 <- summary_disc_sig$beta[i]
    beta_2 <- summary_rep_sig$beta[i]
    sigma_com <- summary_com_sig$se[i]
    sigma_1 <- summary_disc_sig$se[i]
    sigma_2 <- summary_rep_sig$se[i]

    g <- function(x){return(stats::pnorm((x-c_1*sigma_1)/(sigma_1*sigma_com/sigma_2))+stats::pnorm((-x-c_1*sigma_1)/(sigma_1*sigma_com/sigma_2)))}
    f <- function(beta){return((((stats::dnorm((beta_com-beta)/sigma_com))/sigma_com)*g(beta_com))/(stats::pnorm((beta/sigma_1)-c_1)+stats::pnorm(-c_1-(beta/sigma_1))))}

    beta_MLE <- stats::optimize(f, c(sign(beta_com)*max(abs(beta_1),abs(beta_2)),-1*sign(beta_com)*max(abs(beta_1),abs(beta_2))), maximum=TRUE)$maximum
    K <- (sigma_com^2)/((sigma_com^2)+(beta_com-beta_MLE)^2)
    beta_MSE <- K*beta_com + (1-K)*beta_MLE
    beta_adj <- data.frame(beta_MLE,beta_MSE)
    return(beta_adj)
  }

  beta_MLE <- c(rep(0,length(summary_com_sig$beta)))
  beta_MSE <- c(rep(0,length(summary_com_sig$beta)))
  for (b in 1:length(beta_MLE)){
    beta_MLE[b] <- beta_adj(b)$beta_MLE
    beta_MSE[b] <- beta_adj(b)$beta_MSE
    }
  betas <- data.frame(rsid = summary_disc_sig$rsid, beta_disc = summary_disc_sig$beta, se_disc  = summary_disc_sig$se, beta_rep = summary_rep_sig$beta, se_rep = summary_rep_sig$se, beta_com=summary_com_sig$beta,beta_MLE=beta_MLE, beta_MSE=beta_MSE)

  ## reorder based on significance in first set up!
  out <- dplyr::arrange(betas, dplyr::desc(abs(betas$beta_disc/betas$se_disc)))
  return(out)
}

