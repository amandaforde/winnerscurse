#' MSE minimization method for use with discovery and replication GWASs
#'
#' Approach to combine association estimates from discovery and replication data
#' sets
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
#' @param spline Logical parameter, default is spline=FALSE, spline=TRUE uses a
#'   cubic smoothing spline to predict value of B=summary_disc$beta -
#'   summary_rep$beta, spline=FALSE just uses this raw value for B
#'
#' @return Data frame with only SNPs deemed significance in discovery GWAS,
#'   contains their beta and se values from both discovery and replication GWAS
#'   as well as the combined estimate, beta_joint
#' @references Ferguson, J., Alvarez-Iglesias, A., Newell, J., Hinde, J., &
#'   O'Donnell, M. (2017).  Joint incorporation of randomised and observational
#'   evidence in estimating treatment effects. \emph{Statistical Methods in
#'   Medical Research}, \strong{28(1)}, 235\eqn{-}247.
#'   \url{https://doi.org/10.1177/0962280217720854}
#' @export
#'
MSE_minimizer <- function(summary_disc, summary_rep, alpha, spline=FALSE){

  c <- stats::qnorm(1-(alpha)/2)
  summary_disc_sig <- summary_disc[abs(summary_disc$beta/summary_disc$se) > c,]
  summary_rep_sig <- summary_rep[abs(summary_disc$beta/summary_disc$se) > c,]

  B <- summary_disc_sig$beta - summary_rep_sig$beta

  if (spline == FALSE){
    w <- ((summary_rep_sig$se^2)^-1)/(((summary_rep_sig$se^2)^-1) + ((summary_rep_sig$se^2 + B^2)^-1))
  }else{
    z <- summary_disc_sig$beta/summary_disc_sig$se
    B1 <- stats::predict(stats::smooth.spline(x=z,y=B)$fit, z)$y
    w <- ((summary_rep_sig$se^2)^-1)/(((summary_rep_sig$se^2)^-1) + ((summary_rep_sig$se^2 + B1^2)^-1))
  }

  beta_joint <- w*summary_rep_sig$beta + (1-w)*summary_disc_sig$beta
  betas <- data.frame(rsid = summary_disc_sig$rsid, beta_disc = summary_disc_sig$beta, se_disc  = summary_disc_sig$se, beta_rep = summary_rep_sig$beta, se_rep = summary_rep_sig$se, beta_joint=beta_joint)

  ## reorder based on significance in first set up!
  out <- dplyr::arrange(betas, dplyr::desc(abs(betas$beta_disc/betas$se_disc)))
  return(out)
}
