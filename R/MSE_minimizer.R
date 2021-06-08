#' MSE minimization method for use with discovery and replication GWASs
#'
#' \code{MSE_minimizer} is a function which implements an approach that combines
#' the association estimates obtained from discovery and replication GWASs to
#' form a new combined estimate for each SNP. The method used by this function
#' is inspired by that detailed in
#' \href{https://journals.sagepub.com/doi/full/10.1177/0962280217720854}{Ferguson
#' \emph{et al.} (2017)}.
#'
#' @param summary_disc A data frame containing summary statistics from the
#'   \emph{discovery} GWAS. It must have three columns with column names
#'   \code{rsid}, \code{beta} and \code{se}, respectively, and all columns must
#'   contain numerical values. Each row must correspond to a unique SNP,
#'   identified by the numerical value \code{rsid}. The function requires that
#'   there must be at least 5 SNPs as any less will result in issues upon usage
#'   of the smoothing spline.
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
#' @param spline A logical value which determines whether or not a cubic
#'   smoothing spline is to be used. When \code{spline=FALSE}, the value for
#'   \eqn{B} in the formula detailed in the aforementioned paper is merely
#'   calculated as \code{B=summary_disc$beta - summary_rep$beta} for each SNP.
#'   Alternatively, \code{spline=TRUE} applies a cubic smoothing spline to
#'   predict values for \eqn{B} when \code{B=summary_disc$beta -
#'   summary_rep$beta} is regressed on
#'   \code{z=summary_disc$beta/summary_disc$se}, and it is these predicted
#'   values that are then used for \eqn{B}.

#'
#' @return A data frame with summary statistics and adjusted association
#'   estimate of only those SNPs which have been deemed significant in the
#'   discovery GWAS according to the specified threshold, \code{alpha}, i.e.
#'   SNPs with \eqn{p}-values less than \code{alpha}. The inputted summary data
#'   occupies the first five columns, in which the columns \code{beta_disc} and
#'   \code{se_disc} contain the statistics from the discovery GWAS and columns
#'   \code{beta_rep} and \code{se_rep} hold the replication GWAS statistics. The
#'   new combination estimate for each SNPis contained in the final column,
#'   namely \code{beta_joint}. The SNPs are contained in this data frame
#'   according to their significance, with the most significant SNP, i.e. the
#'   SNP with the largest absolute \eqn{z}-statistic, now located in the first
#'   row of the data frame. If no SNPs are detected as significant in the
#'   discovery GWAS, \code{MSE_minimizer} merely returns a data frame which
#'   combines the two inputted data sets.
#'
#' @references Ferguson, J., Alvarez-Iglesias, A., Newell, J., Hinde, J., &
#'   O'Donnell, M. (2017).  Joint incorporation of randomised and observational
#'   evidence in estimating treatment effects. \emph{Statistical Methods in
#'   Medical Research}, \strong{28(1)}, 235\eqn{-}247.
#'   \url{https://doi.org/10.1177/0962280217720854}
#'
#' @seealso
#' \url{https://amandaforde.github.io/winnerscurse/articles/discovery_replication.html}
#' for illustration of the use of \code{MSE_minimizer} with toy data sets and
#' further information regarding computation of the combined SNP-trait
#' association estimates for significant SNPs.
#' @export
#'
MSE_minimizer <- function(summary_disc, summary_rep, alpha=5e-8, spline=TRUE){

  stopifnot(all(c("rsid", "beta","se") %in% names(summary_disc)))
  stopifnot(!all(is.na(summary_disc$rsid)) && !all(is.na(summary_disc$beta)) && !all(is.na(summary_disc$se)))
  stopifnot(is.numeric(summary_disc$rsid) && is.numeric(summary_disc$rsid) && is.numeric(summary_disc$rsid))
  stopifnot(!any(duplicated(summary_disc$rsid)))

  stopifnot(all(c("rsid", "beta","se") %in% names(summary_rep)))
  stopifnot(!all(is.na(summary_rep$rsid)) && !all(is.na(summary_rep$beta)) && !all(is.na(summary_rep$se)))
  stopifnot(is.numeric(summary_rep$rsid) && is.numeric(summary_rep$rsid) && is.numeric(summary_rep$rsid))
  stopifnot(!any(duplicated(summary_rep$rsid)))

  stopifnot(summary_disc$rsid == summary_rep$rsid)
  stopifnot(nrow(summary_disc) > 5)

  c <- stats::qnorm(1-(alpha)/2)

  if(sum(abs(summary_disc$beta/summary_disc$se) > c) == 0){
    summary_data <- cbind(summary_disc[1:3],summary_rep[2:3])
    names(summary_data)[2] <- "beta_disc"
    names(summary_data)[3] <- "se_disc"
    names(summary_data)[4] <- "beta_rep"
    names(summary_data)[5] <- "se_rep"
    return(summary_data)
  }

  B <- summary_disc$beta - summary_rep$beta

  if (spline == FALSE){
    w <- ((summary_rep$se^2)^-1)/(((summary_rep$se^2)^-1) + ((summary_disc$se^2 + B^2)^-1))
  }else{
    z <- summary_disc$beta/summary_disc$se
    B1 <- stats::predict(stats::smooth.spline(x=z,y=B)$fit, z)$y
    w <- ((summary_rep$se^2)^-1)/(((summary_rep$se^2)^-1) + ((summary_disc$se^2 + B1^2)^-1))
  }

  beta_joint <- w*summary_rep$beta + (1-w)*summary_disc$beta

  summary_disc_sig <- summary_disc[abs(summary_disc$beta/summary_disc$se) > c,]
  summary_rep_sig <- summary_rep[abs(summary_disc$beta/summary_disc$se) > c,]
  beta_joint <- beta_joint[abs(summary_disc$beta/summary_disc$se) > c]

  betas <- data.frame(rsid = summary_disc_sig$rsid, beta_disc = summary_disc_sig$beta, se_disc  = summary_disc_sig$se, beta_rep = summary_rep_sig$beta, se_rep = summary_rep_sig$se, beta_joint=beta_joint)

  ## reorder based on significance in first set up!
  out <- dplyr::arrange(betas, dplyr::desc(abs(betas$beta_disc/betas$se_disc)))
  return(out)
}
