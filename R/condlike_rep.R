#' Conditional likelihood method for use with discovery and replication GWASs
#'
#'  Winner's Curse adjustment using discovery and replication summary statistics
#' - Zhong and Prentice (2008)
#'
#' @param summary_disc A data frame containing summary statistics from the
#'   discovery GWAS. It must have three columns with column names \code{rsid},
#'   \code{beta} and \code{se}, respectively, and all columns must contain
#'   numerical values.
#' @param summary_rep A data frame containing summary statistics from the
#'   replication GWAS. It must have three columns with column names \code{rsid},
#'   \code{beta} and \code{se}, respectively, and all columns must contain
#'   numerical values.
#' @param alpha A numerical value which specifies the desired genome-wide
#'   significance threshold for the discovery GWAS.
#' @param conf_interval A logical value.
#' @param conf_level A numerical value.
#'
#' @return Data frame with only SNPs deemed significance in discovery GWAS,
#'   contains their beta and se values from both discovery and replication GWAS
#'   as well as the combined unadjusted estimate, the conditional likelihood MLE
#'   and the weighted average for each SNP
#'
#' @references Zhong, H., & Prentice, R. L. (2008). Bias-reduced estimators and
#'   confidence intervals for odds ratios in genome-wide association studies.
#'   \emph{Biostatistics (Oxford, England)}, \strong{9(4)}, 621\eqn{-}634.
#'   \url{https://doi.org/10.1093/biostatistics/kxn001}
#' @export
#'
condlike_rep <- function(summary_disc,summary_rep,alpha, conf_interval=FALSE, conf_level=0.95){
  ### only selection on snps from discovery data set, just require definition for c or equivalently, alpha

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
    logf <- function(beta){return(log(f(beta)))}

    beta_MLE <- stats::optimize(f, c(sign(beta_com)*max(abs(beta_1),abs(beta_2)),-1*sign(beta_com)*max(abs(beta_1),abs(beta_2))), maximum=TRUE)$maximum

    fun_CI <- function(beta){return(logf(beta) - logf(beta_MLE) + (qchisq(conf_level, df=1))/2)}
    beta_MLE_lower <- uniroot(fun_CI, interval=c(beta_MLE+round(min(summary_disc$beta),1),beta_MLE))$root
    beta_MLE_upper <- uniroot(fun_CI, interval=c(beta_MLE,beta_MLE+round(max(summary_disc$beta),1)))$root
    beta_adj <- data.frame(beta_MLE,beta_MLE_lower, beta_MLE_upper)
    return(beta_adj)
  }

  beta_MLE <- c(rep(0,length(summary_com_sig$beta)))
  for (b in 1:length(beta_MLE)){
    beta_MLE[b] <- beta_adj(b)$beta_MLE
  }

  K <- (summary_com_sig$se^2)/((summary_com_sig$se^2)+(summary_com_sig$beta-beta_MLE)^2)
  beta_MSE <- K*summary_com_sig$beta + (1-K)*beta_MLE

  betas <- data.frame(rsid = summary_disc_sig$rsid, beta_disc = summary_disc_sig$beta, se_disc  = summary_disc_sig$se, beta_rep = summary_rep_sig$beta, se_rep = summary_rep_sig$se, beta_com=summary_com_sig$beta,beta_MLE=beta_MLE, beta_MSE=beta_MSE)

  ## reorder based on significance in first set up!
  out <- dplyr::arrange(betas, dplyr::desc(abs(betas$beta_disc/betas$se_disc)))

  if (conf_interval == FALSE){return(out)}else{
    beta_com_lower <- summary_com_sig$beta - qnorm(1-(1-conf_level)/2)*summary_com_sig$se
    beta_com_upper <- summary_com_sig$beta + qnorm(1-(1-conf_level)/2)*summary_com_sig$se

    beta_MLE_lower <- c(rep(0,length(summary_com_sig$beta)))
    beta_MLE_upper <- c(rep(0,length(summary_com_sig$beta)))
    for (b in 1:length(beta_MLE_lower)){
      beta_MLE_lower[b] <- beta_adj(b)$beta_MLE_lower
      beta_MLE_upper[b] <- beta_adj(b)$beta_MLE_upper
    }

    K_lower <- (summary_com_sig$se^2)/((summary_com_sig$se^2)+(beta_com_lower-beta_MLE_lower)^2)
    beta_MSE_lower <- K_lower*beta_com_lower + (1-K_lower)*beta_MLE_lower

    K_upper <- (summary_com_sig$se^2)/((summary_com_sig$se^2)+(beta_com_upper-beta_MLE_upper)^2)
    beta_MSE_upper <- K_upper*beta_com_upper + (1-K_upper)*beta_MLE_upper

    betas <- data.frame(rsid = summary_disc_sig$rsid, beta_disc = summary_disc_sig$beta, se_disc  = summary_disc_sig$se, beta_rep = summary_rep_sig$beta, se_rep = summary_rep_sig$se, beta_com=summary_com_sig$beta, beta_com_lower=beta_com_lower, beta_com_upper=beta_com_upper, beta_MLE = beta_MLE, beta_MLE_lower=beta_MLE_lower, beta_MLE_upper = beta_MLE_upper, beta_MSE, beta_MSE_lower,beta_MSE_upper)
    ## reorder based on significance in first set up!
    out <- dplyr::arrange(betas, dplyr::desc(abs(betas$beta_disc/betas$se_disc)))
    return(out)


    }
}

