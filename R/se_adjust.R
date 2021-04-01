#' Standard errors of adjusted discovery GWAS estimates via parametric bootstrap
#'
#' \code{se_adjust} is a function which allows the user to obtain approximate
#' standard errors of adjusted association estimates, by means of parametric
#' bootstrapping. Standard errors can be evaluated for estimates which have been
#' corrected with the Empirical Bayes method, FDR Inverse Quantile
#' Transformation method or the bootstrap method. Note that in comparison to the
#' other functions in this package, this function can be computationally
#' intensive and take a several minutes to run, depending on the size of the
#' data set and the number of bootstraps chosen.
#'
#' @param summary_data A data frame containing summary statistics from the
#'   discovery GWAS. It must have three columns with column names \code{rsid},
#'   \code{beta} and \code{se}, respectively, and all columns must contain
#'   numerical values.
#' @param method A string specifying the function to be implemented on each of
#'   the bootstrap samples. It should take the form \code{"BR_ss"},
#'   \code{"empirical_bayes"} or \code{"FDR_IQT"}.
#' @param n_boot A numerical value which determines the number of bootstrap
#'   repetitions to be used. The default value is \code{100}.
#'
#' @return A data frame which combines the output of the chosen method with an
#'   additional column, namely \code{adj_se}. This column provides the standard
#'   errors of the adjusted association estimates for each SNP.
#' @seealso \code{\link{empirical_bayes}}, \code{\link{BR_ss}} and
#'   \code{\link{FDR_IQT}} for details on operation of these methods with
#'   summary statistics from discovery GWAS.
#'
#'   \url{https://amandaforde.github.io/winnerscurse/articles/standard_errors_confidence_intervals.html}
#'    for illustration of the use of \code{se_adjust} with a toy data set and
#'   further information regarding the manner in which the standard errors are
#'   computed.
#' @export
#'
se_adjust <- function(summary_data, method, n_boot = 100){
  N <- nrow(summary_data)
  beta_boot <- matrix(stats::rnorm(n_boot*N, mean = rep(summary_data$beta,n_boot), sd = rep(summary_data$se,n_boot)), nrow=N, ncol=n_boot, byrow=FALSE)
  beta_adj_boot <- matrix(nrow=N, ncol=n_boot)

  for (i in 1:n_boot){
    stats_boot <- data.frame(rsid = summary_data$rsid, beta = beta_boot[,i], se = summary_data$se)
    if (method == "BR_ss"){out <- BR_ss(stats_boot)}
    if (method == "FDR_IQT") {out <- FDR_IQT(stats_boot)}
    if (method == "empirical_bayes") {out <- empirical_bayes(stats_boot)}
    out <- dplyr::arrange(out,out$rsid)
    beta_adj_boot[,i] <- out[,4]
  }

  adj_se <- apply(beta_adj_boot,1,stats::sd)
  data_se <- cbind(summary_data, adj_se)
  data_se <- dplyr::arrange(data_se, dplyr::desc(abs(data_se$beta/data_se$se)))
  if (method == "BR_ss"){summary_data <- BR_ss(summary_data)}
  if (method == "FDR_IQT"){summary_data <- FDR_IQT(summary_data)}
  if (method == "empirical_bayes"){summary_data <- empirical_bayes(summary_data)}
  summary_data <- cbind(summary_data, data_se[4])
  return(summary_data)
}


