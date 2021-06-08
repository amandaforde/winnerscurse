#'Bootstrap method for use with discovery GWAS
#'
#'\code{BR_ss} is a function which aims to use summary statistics to alleviate
#'Winner's Curse bias in SNP-trait association estimates, obtained from a
#'discovery GWAS. The function implements an adaptation of a bootstrap
#'resampling method known as BR-squared, detailed in
#'\href{https://onlinelibrary.wiley.com/doi/abs/10.1002/sim.4228}{Faye \emph{et
#'al.} (2011)}.
#'
#'@param summary_data A data frame containing summary statistics from the
#'  discovery GWAS. It must have three columns with column names \code{rsid},
#'  \code{beta} and \code{se}, respectively, and all columns must contain
#'  numerical values. Each row must correspond to a unique SNP, identified by
#'  the numerical value \code{rsid}. The function requires that there must be at
#'  least 5 SNPs as any less will result in issues upon usage of the
#'  smoothing spline.
#'@param seed_opt A logical value which allows the user to choose if they wish
#'  to set a seed, in order to ensure reproducibility of adjusted estimates.
#'  Small differences can occur between iterations of the function with the same
#'  data set due to the use of parametric bootstrapping. The default setting is
#'  \code{seed_opt=FALSE}.
#'@param seed A numerical value which specifies the seed used if
#'  \code{seed_opt=TRUE}. The default setting is the arbitrary value of
#'  \code{1998}.
#'
#'@return A data frame with the inputted summary data occupying the first three
#'  columns. The new adjusted association estimates for each SNP are returned in
#'  the fourth column, namely \code{beta_BR_ss}. The SNPs are contained in this
#'  data frame according to their significance, with the most significant SNP,
#'  i.e. the SNP with the largest absolute \eqn{z}-statistic, now located in the
#'  first row of the data frame.
#'@references Faye, L. L., Sun, L., Dimitromanolakis, A., & Bull, S. B. (2011).
#'  A flexible genome-wide bootstrap method that accounts for ranking and
#'  threshold-selection bias in GWAS interpretation and replication study
#'  design. \emph{Statistics in Medicine}, \strong{30(15)}, 1898\eqn{-}1912.
#'  \url{https://doi.org/10.1002/sim.4228}
#'
#'@seealso
#'\url{https://amandaforde.github.io/winnerscurse/articles/winners_curse_methods.html}
#'for illustration of the use of \code{BR_ss} with a toy data set and further
#'information regarding the computation of the adjusted SNP-trait association
#'estimates.
#'
#'
#'@export

BR_ss <- function(summary_data,seed_opt = FALSE,seed=1998){

    stopifnot(all(c("rsid", "beta","se") %in% names(summary_data)))
    stopifnot(!all(is.na(summary_data$rsid)) && !all(is.na(summary_data$beta)) && !all(is.na(summary_data$se)))
    stopifnot(nrow(summary_data) > 5)
    stopifnot(is.numeric(summary_data$rsid) && is.numeric(summary_data$rsid) && is.numeric(summary_data$rsid))
    stopifnot(!any(duplicated(summary_data$rsid)))

    summary_data <- dplyr::arrange(summary_data, dplyr::desc((summary_data$beta/summary_data$se)))

    N <- nrow(summary_data)
    if(seed_opt==TRUE){set.seed(seed)}
    beta_boot <- matrix(stats::rnorm(1*N, mean = rep(summary_data$beta,1), sd = rep(summary_data$se,1)), nrow=N, ncol=1, byrow=FALSE)

    beta_mat <- matrix(rep(summary_data$beta,1), nrow=N, ncol=1, byrow=FALSE)
    se_mat <- matrix(rep(summary_data$se,1), nrow=N, ncol=1, byrow=FALSE)

    beta_oob <- beta_mat

    ordering <- apply(beta_boot/se_mat, 2,order,decreasing=TRUE)
    bias_correct <- matrix(nrow=N, ncol=1)


    bias_correct[,1] <- (beta_boot[ordering[,1],1] - beta_oob[ordering[,1],1])/summary_data$se[ordering[,1]]


    z <- summary_data$beta/summary_data$se

    bias_correct <- stats::predict(stats::smooth.spline(z,bias_correct)$fit, z)$y

    beta_BR_ss <- summary_data$beta - summary_data$se*bias_correct[rank(-1*summary_data$beta/summary_data$se)]
    beta_BR_ss[sign(beta_BR_ss) != sign(summary_data$beta)] <- 0

    summary_data <- cbind(summary_data, beta_BR_ss)

    for (i in 1:N){
      if(abs(summary_data$beta[i]) < abs(summary_data$beta_BR_ss[i])){summary_data$beta_BR_ss[i] <- summary_data$beta[i]}
    }

    summary_data <- dplyr::arrange(summary_data,dplyr::desc(abs(summary_data$beta/summary_data$se)))

    return(summary_data)
  }
