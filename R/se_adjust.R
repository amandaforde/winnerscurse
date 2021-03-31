#' Standard errors of adjusted discovery GWAS estimates via parametric bootstrap
#'
#' @param summary_data A data frame containing summary statistics from the
#'   discovery GWAS. It must have three columns with column names \code{rsid},
#'   \code{beta} and \code{se}, respectively, and all columns must contain
#'   numerical values.
#' @param method A string which specifies the function to be implemented
#'   \code{"BR_ss"}, \code{"empirical_bayes"} \code{"FDR_IQT"}
#' @param n_boot A numerical value - number of bootstrap repetitions used -
#'   defaults to 100
#'
#' @return Data frame with summary data together with adjusted estimates and
#'   standard error: adj_se
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


