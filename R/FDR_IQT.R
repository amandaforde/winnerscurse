#' FDR Inverse Quantile Transformation method to correct for Winner’s Curse - Bigdeli et al. (2016)
#'
#' @param summary_data Data frame containing summary data, three columns: rsid, beta, se
#' @param min_pval avoids zero p-values (in R) which gives trouble with qnorm(), also abs(z) > 8 is generally not biased
#'
#' @return Data frame with summary data together with adjusted estimates
#' @export
#'
#'
FDR_IQT <- function(summary_data, min_pval=1e-15)

{
  z <- summary_data$beta/summary_data$se
  p_val <- 2*(1-stats::pnorm(abs(z)))
  p_val[p_val < min_pval] <- min_pval

  adj_p <- stats::p.adjust(p_val, method="fdr")
  adj_z <- stats::qnorm(1-(adj_p/2))
  adj_z[abs(z) > stats::qnorm(1-(min_pval/2))] <- z[abs(z) > stats::qnorm(1-(min_pval/2))]

  beta_FIQT <- sign(summary_data$beta)*adj_z*summary_data$se
  summary_data <- cbind(summary_data,beta_FIQT)

  summary_data <- dplyr::arrange(summary_data,dplyr::desc(abs(summary_data$beta/summary_data$se)))

  return(summary_data)

}
