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
#'  presumption that in general, estimates of SNPs with \eqn{z > 8} are not
#'  biased. The default value is \code{min_pval = 1e-15}.
#'@param method A string which allows the user to choose if they wish to execute
#'  the original definition of the method or an alternative form which adjusts
#'  the effect sizes of all SNPs, not only those that have \eqn{p}-values
#'  greater than \code{1e-15}. The default setting is \code{method="original"}
#'  while the other form of the method can be executed using
#'  \code{method="alt"}.
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
#'  \url{https://doi.org/10.1093/bioinformatics/btw303}
#'
#'@seealso
#'\url{https://amandaforde.github.io/winnerscurse/articles/winners_curse_methods.html}
#'for illustration of the use of \code{FDR_IQT} with a toy data set and further
#'information regarding the computation of the adjusted SNP-trait association
#'estimates.
#'@export
#'
#'
FDR_IQT <- function(summary_data, min_pval=1e-15, method="original")

{

  stopifnot(all(c("rsid", "beta","se") %in% names(summary_data)))
  stopifnot(!all(is.na(summary_data$rsid)) && !all(is.na(summary_data$beta)) && !all(is.na(summary_data$se)))
  stopifnot(is.numeric(summary_data$beta) && is.numeric(summary_data$se))
  stopifnot(!any(duplicated(summary_data$rsid)))


  z <- summary_data$beta/summary_data$se
  p_val <- 2*(1-stats::pnorm(abs(z)))

  if(method=="original"){
    p_val[p_val < min_pval] <- min_pval
    adj_p <- stats::p.adjust(p_val, method="fdr")
    adj_z <- stats::qnorm(1-(adj_p/2))
    adj_z[abs(z) > stats::qnorm(1-(min_pval/2))] <- abs(z[abs(z) > stats::qnorm(1-(min_pval/2))])
    beta_FIQT <- sign(summary_data$beta)*adj_z*summary_data$se
  }


  if(method=="alt"){
    p_val[p_val < min_pval] <- 2*(exp(-abs(z[p_val < min_pval])^2/2)/(abs(z[p_val < min_pval])*sqrt(2*pi)))
    adj_p <- stats::p.adjust(p_val, method="fdr")
    adj_z <- stats::qnorm(1-(adj_p/2))
    adjusted <- data.frame(adj_p,adj_z)
    zz <- seq(from = 8, to = floor(max(z)) + 1, by = 0.0001)
    p_vals <- 2*(exp(-abs(zz)^2/2)/(abs(zz)*sqrt(2*pi)))
    for (i in 1:length(p_val[p_val < min_pval])){
      if(adjusted$adj_z[i] == Inf){
        adjusted$adj_z[i] <- zz[which.min(abs(p_vals - c(rep(adjusted$adj_p[i],length(p_vals)))))]
      }
    }
    beta_FIQT <- sign(summary_data$beta)*adjusted$adj_z*summary_data$se
  }

  summary_data <- cbind(summary_data,beta_FIQT)
  summary_data <- dplyr::arrange(summary_data,dplyr::desc(abs(summary_data$beta/summary_data$se)))

  return(summary_data)

}


