#'FDR IQT method for use with discovery GWAS
#'
#'\code{FDR_IQT} is a function which uses summary statistics to reduce Winner's
#'Curse bias in SNP-trait association estimates, obtained from a discovery GWAS.
#'The function implements the FDR Inverse Quantile Transformation method
#'described in
#'\href{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5013908/}{Bigdeli \emph{et
#'al.} (2016)}, which was established for this purpose.
#'@param summary_data A data frame containing summary statistics from the
#'  discovery GWAS. It must have three columns with column names \code{rsid},
#'  \code{beta} and \code{se}, respectively, and columns \code{beta} and
#'  \code{se} must contain numerical values. Each row must correspond to a
#'  unique SNP, identified by the numerical value \code{rsid}.
#'@param min_pval A numerical value whose purpose is to avoid zero
#'  \eqn{p}-values as this introduces issues when \code{qnorm()} is applied. Any
#'  SNP for which its computed \eqn{p}-value is found to be less than
#'  \code{min_pval} is merely re-assigned \code{min_pval} as its \eqn{p}-value
#'  and the method proceeds. By definition, the method makes no adjustment to
#'  the association estimate of a SNP for which this has occurred with the
#'  presumption that in general, estimates of SNPs with \eqn{z > 37} are not
#'  biased. The default value is \code{min_pval = 1e-300}.
#'
#'@return A data frame with the inputted summary data occupying the first three
#'  columns. The new adjusted association estimates for each SNP are returned in
#'  the fourth column, namely \code{beta_FIQT}. The SNPs are contained in this
#'  data frame according to their significance, with the most significant SNP,
#'  i.e. the SNP with the largest absolute \eqn{z}-statistic, now located in the
#'  first row of the data frame.
#'
#'
#'@references Bigdeli, T. B., Lee, D., Webb, B. T., Riley, B. P., Vladimirov, V.
#'  I., Fanous, A. H., Kendler, K. S., & Bacanu, S. A. (2016). A simple yet
#'  accurate correction for winner's curse can predict signals discovered in
#'  much larger genome scans. \emph{Bioinformatics (Oxford, England)},
#'  \strong{32(17)}, 2598\eqn{-}2603.
#'  \doi{10.1093/bioinformatics/btw303}
#'
#'@seealso
#'\url{https://amandaforde.github.io/winnerscurse/articles/winners_curse_methods.html}
#'for illustration of the use of \code{FDR_IQT} with a toy data set and further
#'information regarding the computation of the adjusted SNP-trait association
#'estimates.
#'@export
#'
#'
FDR_IQT <- function(summary_data, min_pval=1e-300){

  stopifnot(all(c("rsid", "beta","se") %in% names(summary_data)))
  stopifnot(!all(is.na(summary_data$rsid)) && !all(is.na(summary_data$beta)) && !all(is.na(summary_data$se)))
  stopifnot(is.numeric(summary_data$beta) && is.numeric(summary_data$se))
  stopifnot(!any(duplicated(summary_data$rsid)))

  z <- summary_data$beta/summary_data$se
  p_val <- 2*(stats::pnorm(abs(z), lower.tail=FALSE))

  p_val[p_val < min_pval] <- min_pval
  adj_p <- stats::p.adjust(p_val, method="fdr")
  adj_z <- stats::qnorm((adj_p/2), lower.tail=FALSE)
  adj_z[abs(z) > stats::qnorm((min_pval/2), lower.tail=FALSE)] <- abs(z[abs(z) > stats::qnorm((min_pval/2), lower.tail=FALSE)])
  beta_FIQT <- sign(summary_data$beta)*adj_z*summary_data$se


  summary_data <- cbind(summary_data,beta_FIQT)
  summary_data <- dplyr::arrange(summary_data,dplyr::desc(abs(summary_data$beta/summary_data$se)))

  return(summary_data)
}


