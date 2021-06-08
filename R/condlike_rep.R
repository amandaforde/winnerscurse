#' Conditional likelihood method for use with discovery and replication GWASs
#'
#' \code{condlike_rep} is a function which attempts to produce less biased
#' SNP-trait association estimates for SNPs deemed significant in the discovery
#' GWAS, using summary statistics from both discovery and replication GWASs. The
#' function computes three new association estimates for each SNP in a manner
#' based closely on the method described in
#' \href{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2536726/}{Zhong and
#' Prentice (2008)}. It also returns confidence intervals for each new
#' association estimate, if desired by the user.
#'
#' @param summary_disc A data frame containing summary statistics from the
#'   \emph{discovery} GWAS. It must have three columns with column names
#'   \code{rsid}, \code{beta} and \code{se}, respectively, and all columns must
#'   contain numerical values. Each row must correspond to a unique SNP,
#'   identified by the numerical value \code{rsid}.
#' @param summary_rep A data frame containing summary statistics from the
#'   \emph{replication} GWAS. It must have three columns with column names
#'   \code{rsid}, \code{beta} and \code{se}, respectively, and all columns must
#'   contain numerical values. Each row must correspond to a unique SNP,
#'   identified by the numerical value \code{rsid}. SNPs must be ordered in the
#'   exact same manner as those in \code{summary_disc}, i.e.
#'   \code{summary_rep$rsid} must be equivalent to \code{summary_disc$rsid}.
#' @param alpha A numerical value which specifies the desired genome-wide
#'   significance threshold for the discovery GWAS. The default is given as
#'   \code{5e-8}.
#' @param conf_interval A logical value which determines whether or not
#'   confidence intervals for each form of adjusted association estimate is also
#'   to be computed and outputted. The default is \code{conf_interval=FALSE}.
#' @param conf_level A numerical value between 0 and 1 which specifies the
#'   confidence interval to be computed. The default setting is \code{0.95}
#'   which results in the calculation of a 95\% confidence interval for the
#'   adjusted association estimate for each SNP.
#'
#' @return A data frame with summary statistics and adjusted association
#'   estimates of only those SNPs which have been deemed significant in the
#'   discovery GWAS according to the specified threshold, \code{alpha}, i.e.
#'   SNPs with \eqn{p}-values less than \code{alpha}. The inputted summary data
#'   occupies the first five columns, in which the columns \code{beta_disc} and
#'   \code{se_disc} contain the statistics from the discovery GWAS and columns
#'   \code{beta_rep} and \code{se_rep} hold the replication GWAS statistics. For
#'   the default setting of \code{conf_interval=FALSE}, the new adjusted
#'   association estimates for each SNP, as defined in the aforementioned paper,
#'   are contained in the next three columns, namely \code{beta_com},
#'   \code{beta_MLE} and \code{beta_MSE}. For the case when
#'   \code{conf_interval=TRUE}, the lower and upper boundaries of each
#'   confidence interval for each form of adjusted estimate are included in the
#'   data frame as well as the adjusted estimates for each SNP. The SNPs are
#'   contained in this data frame according to their significance, with the most
#'   significant SNP, i.e. the SNP with the largest absolute \eqn{z}-statistic,
#'   now located in the first row of the data frame. If no SNPs are detected as
#'   significant in the discovery GWAS, \code{condlike_rep} merely returns a
#'   data frame which combines the two inputted data sets.
#'
#' @references Zhong, H., & Prentice, R. L. (2008). Bias-reduced estimators and
#'   confidence intervals for odds ratios in genome-wide association studies.
#'   \emph{Biostatistics (Oxford, England)}, \strong{9(4)}, 621\eqn{-}634.
#'   \url{https://doi.org/10.1093/biostatistics/kxn001}
#'
#' @seealso
#'
#' \url{https://amandaforde.github.io/winnerscurse/articles/discovery_replication.html}
#' for illustration of the use of \code{condlike_rep} with toy data sets and
#' further information regarding computation of the adjusted SNP-trait
#' association estimates and their corresponding confidence intervals for
#' significant SNPs.
#' @export
#'
condlike_rep <- function(summary_disc,summary_rep,alpha=5e-8, conf_interval=FALSE, conf_level=0.95){
  ### only selection on snps from discovery data set, just require definition for c or equivalently, alpha

  stopifnot(all(c("rsid", "beta","se") %in% names(summary_disc)))
  stopifnot(!all(is.na(summary_disc$rsid)) && !all(is.na(summary_disc$beta)) && !all(is.na(summary_disc$se)))
  stopifnot(is.numeric(summary_disc$rsid) && is.numeric(summary_disc$rsid) && is.numeric(summary_disc$rsid))
  stopifnot(!any(duplicated(summary_disc$rsid)))

  stopifnot(all(c("rsid", "beta","se") %in% names(summary_rep)))
  stopifnot(!all(is.na(summary_rep$rsid)) && !all(is.na(summary_rep$beta)) && !all(is.na(summary_rep$se)))
  stopifnot(is.numeric(summary_rep$rsid) && is.numeric(summary_rep$rsid) && is.numeric(summary_rep$rsid))
  stopifnot(!any(duplicated(summary_rep$rsid)))

  stopifnot(summary_disc$rsid == summary_rep$rsid)

  se_com <- sqrt((((summary_disc$se)^2)*((summary_rep$se)^2))/(((summary_disc$se)^2) + ((summary_rep$se)^2)))
  beta_com <- ((((summary_rep$se)^2)*(summary_disc$beta))+(((summary_disc$se)^2)*(summary_rep$beta)))/(((summary_disc$se)^2) + ((summary_rep$se)^2))
  summary_com <- data.frame(rsid=summary_disc$rsid,beta=beta_com,se=se_com)
  c_1 <- stats::qnorm(1-(alpha)/2)

  if(sum(abs(summary_disc$beta/summary_disc$se) > c_1) == 0){
    summary_data <- cbind(summary_disc[1:3],summary_rep[2:3])
    names(summary_data)[2] <- "beta_disc"
    names(summary_data)[3] <- "se_disc"
    names(summary_data)[4] <- "beta_rep"
    names(summary_data)[5] <- "se_rep"
    return(summary_data)
  }

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

    if (conf_interval == FALSE){beta_adj <- data.frame(beta_MLE)}else{
        fun_CI <- function(beta){return(logf(beta) - logf(beta_MLE) + (stats::qchisq(conf_level, df=1))/2)}
        beta_MLE_lower <- stats::uniroot(fun_CI, interval=c(beta_MLE+round(min(summary_disc$beta),1),beta_MLE))$root
        beta_MLE_upper <- stats::uniroot(fun_CI, interval=c(beta_MLE,beta_MLE+round(max(summary_disc$beta),1)))$root
        beta_adj <- data.frame(beta_MLE,beta_MLE_lower, beta_MLE_upper)
      }
    return(beta_adj)

  }

  beta_MLE <- c(rep(0,length(summary_com_sig$beta)))
  for (b in 1:length(beta_MLE)){
    beta_MLE[b] <- abs(beta_adj(b)$beta_MLE)*sign(summary_com_sig$beta[b])
  }

  K <- (summary_com_sig$se^2)/((summary_com_sig$se^2)+(summary_com_sig$beta-beta_MLE)^2)
  beta_MSE <- K*summary_com_sig$beta + (1-K)*beta_MLE

  betas <- data.frame(rsid = summary_disc_sig$rsid, beta_disc = summary_disc_sig$beta, se_disc  = summary_disc_sig$se, beta_rep = summary_rep_sig$beta, se_rep = summary_rep_sig$se, beta_com=summary_com_sig$beta,beta_MLE=beta_MLE, beta_MSE=beta_MSE)

  ## reorder based on significance in first set up!
  out <- dplyr::arrange(betas, dplyr::desc(abs(betas$beta_disc/betas$se_disc)))

  if (conf_interval == FALSE){return(out)}else{
    beta_com_lower <- summary_com_sig$beta - stats::qnorm(1-(1-conf_level)/2)*summary_com_sig$se
    beta_com_upper <- summary_com_sig$beta + stats::qnorm(1-(1-conf_level)/2)*summary_com_sig$se

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

