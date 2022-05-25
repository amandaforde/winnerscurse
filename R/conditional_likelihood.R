#' Conditional likelihood methods for use with discovery GWAS
#'
#' \code{conditional_likelihood} is a function which uses summary statistics to
#' correct bias created by the Winner's Curse phenomenon in the SNP-trait
#' association estimates, obtained from a discovery GWAS, of SNPs which are
#' considered significant. The function implements the approximate conditional
#' likelihood approach, discussed in
#' \href{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2665019/}{Ghosh \emph{et
#' al.} (2008)}, which suggests three different forms of a less biased
#' association estimate. Note that if the \eqn{z}-statistic of a particular SNP
#' is greater than 100, then merely the original naive estimate will be
#' outputted for the second form of the adjusted estimate, namely
#' \code{beta.cl2}, for that SNP.
#'
#' @param summary_data A data frame containing summary statistics from the
#'   discovery GWAS. It must have three columns with column names \code{rsid},
#'   \code{beta} and \code{se}, respectively, and columns \code{beta} and
#'   \code{se} must contain numerical values. Each row must correspond to a
#'   unique SNP, identified by \code{rsid}.
#' @param alpha A numerical value which specifies the desired genome-wide
#'   significance threshold. The default is given as \code{5e-8}.
#'
#' @return A data frame with summary statistics and adjusted association
#'   estimates of only those SNPs which have been deemed significant according
#'   to the specified threshold, \code{alpha}, i.e. SNPs with \eqn{p}-values
#'   less than \code{alpha}. The inputted summary data occupies the first three
#'   columns. The new adjusted association estimates for each SNP, as defined in
#'   the aforementioned paper, are contained in the next three columns, namely
#'   \code{beta.cl1}, \code{beta.cl2} and \code{beta.cl3}. The SNPs are
#'   contained in this data frame according to their significance, with the most
#'   significant SNP, i.e. the SNP with the largest absolute \eqn{z}-statistic,
#'   now located in the first row of the data frame. However, if no SNPs are
#'   detected as significant in the data set, \code{conditional_likelihood}
#'   returns a warning message: \code{"WARNING: There are no significant SNPs at
#'   this threshold."}
#'
#'
#' @references Ghosh, A., Zou, F., & Wright, F. A. (2008). Estimating odds
#'   ratios in genome scans: an approximate conditional likelihood approach.
#'   \emph{American journal of human genetics}, \strong{82(5)}, 1064\eqn{-}1074.
#'   \url{https://doi.org/10.1016/j.ajhg.2008.03.002}
#'
#'
#' @seealso
#' \url{https://amandaforde.github.io/winnerscurse/articles/winners_curse_methods.html}
#' for illustration of the use of \code{conditional_likelihood} with a toy data
#' set and further information regarding the computation of the adjusted
#' SNP-trait association estimates for significant SNPs.
#' @export
#'
#'
conditional_likelihood <- function(summary_data, alpha=5e-8){

  stopifnot(all(c("rsid", "beta","se") %in% names(summary_data)))
  stopifnot(!all(is.na(summary_data$rsid)) && !all(is.na(summary_data$beta)) && !all(is.na(summary_data$se)))
  stopifnot(is.numeric(summary_data$beta) && is.numeric(summary_data$se))
  stopifnot(!any(duplicated(summary_data$rsid)))


  z <- summary_data$beta/summary_data$se
  p_val <- 2*(1-stats::pnorm(abs(z)))
  summary_data <- cbind(summary_data, z, p_val)


  if (sum(summary_data$p_val<alpha) == 0) {
    warn <- "WARNING: There are no significant SNPs at this threshold."
    return(warn)
  }

  summary_data_sig <- summary_data[summary_data$p_val<alpha,]

  c <- stats::qnorm(1-(alpha)/2)

  beta.cl1 <- c(rep(0,nrow(summary_data_sig)))
  beta.cl2 <- c(rep(0,nrow(summary_data_sig)))
  beta.cl3 <- c(rep(0,nrow(summary_data_sig)))

  for (i in 1:nrow(summary_data_sig)){
    cond.like <- function (mu) { return ((stats::dnorm(summary_data_sig$z[i]-mu))/(stats::pnorm(-c+mu)+stats::pnorm(-c-mu)))}
    mean.cond.like <- function (mu) {return(mu*cond.like(mu))}

    beta.cl1[i] <- (stats::optimize(cond.like, c(0,summary_data_sig$z[i]), maximum=TRUE)$maximum)*summary_data_sig$se[i]

    if(abs(summary_data_sig$z[i]) < 37){
      beta.cl2[i] <- ((stats::integrate(mean.cond.like,-37,37)$value)/(stats::integrate(cond.like,-37,37)$value))*(summary_data_sig$se[i])
    }else{
      if(abs(summary_data_sig$z[i]) > 100){
        beta.cl2[i] <- summary_data_sig$beta[i]
      }else{
        beta.cl2[i] <- ((stats::integrate(mean.cond.like,-100,100)$value)/(stats::integrate(cond.like,-100,100)$value))*(summary_data_sig$se[i])
      }
    }

    beta.cl3[i] <- (beta.cl1[i]+beta.cl2[i])/2
  }

  summary_data_sig <- cbind(summary_data_sig,beta.cl1,beta.cl2,beta.cl3)


  if(sum(abs(beta.cl2) > abs(summary_data_sig$beta + sign(summary_data_sig$beta))) >= 1){
    bad.cl2 <- which((abs(beta.cl2) > abs(summary_data_sig$beta + sign(summary_data_sig$beta))) == TRUE)
    bad.cl2.df <- data.frame(rsid = summary_data_sig$rsid[bad.cl2], z = summary_data_sig$z[bad.cl2], z.cl2 = beta.cl2[bad.cl2]/summary_data_sig$se[bad.cl2], beta = summary_data_sig$beta[bad.cl2], beta.cl2 = beta.cl2[bad.cl2], se = summary_data_sig$se[bad.cl2])
    good.cl2.df <- data.frame(rsid = summary_data_sig$rsid[-(bad.cl2)], z = summary_data_sig$z[-(bad.cl2)], z.cl2 = beta.cl2[-(bad.cl2)]/summary_data_sig$se[-(bad.cl2)], beta = summary_data_sig$beta[-bad.cl2], beta.cl2 = beta.cl2[-(bad.cl2)], se = summary_data_sig$se[-(bad.cl2)])
    for(i in 1:nrow(bad.cl2.df)){
      if(length(good.cl2.df$z[good.cl2.df$z >= bad.cl2.df$z[i]]) == 0 || length(good.cl2.df$z[good.cl2.df$z < bad.cl2.df$z[i]]) == 0){bad.cl2.df$z.cl2[i] <- good.cl2.df$z.cl2[which.min(abs(good.cl2.df$z-bad.cl2.df$z[i]))]
      bad.cl2.df$beta.cl2[i] <- bad.cl2.df$z.cl2[i]*bad.cl2.df$se[i]}else{
        min_greater <- min(good.cl2.df$z[good.cl2.df$z >= bad.cl2.df$z[i]])
        min_smaller <- min(good.cl2.df$z[good.cl2.df$z < bad.cl2.df$z[i]])
        bad.cl2.df$z.cl2[i] <- ((abs(min_greater - bad.cl2.df$z[i]))*(good.cl2.df$z.cl2[good.cl2.df$z == min_smaller]) + (abs(bad.cl2.df$z[i] -     min_smaller))*(good.cl2.df$z.cl2[good.cl2.df$z == min_greater]))/(abs(min_greater - min_smaller))
        bad.cl2.df$beta.cl2[i] <- bad.cl2.df$z.cl2[i]*bad.cl2.df$se[i]
      }

      summary_data_sig$beta.cl2[summary_data_sig$rsid == bad.cl2.df$rsid[i]] <- bad.cl2.df$beta.cl2[i]
      summary_data_sig$beta.cl3[summary_data_sig$rsid == bad.cl2.df$rsid[i]] <- (summary_data_sig$beta.cl1[summary_data_sig$rsid == bad.cl2.df$rsid[i]] +bad.cl2.df$beta.cl2[i])/2
    }
  }

  summary_data_sig <- cbind(summary_data_sig,beta.cl1,beta.cl2,beta.cl3)
  summary_data_sig <- dplyr::arrange(summary_data_sig,dplyr::desc(abs(summary_data_sig$z)))
  summary_data_sig <- cbind(summary_data_sig[,1:3],summary_data_sig[,6:8])
  return(summary_data_sig)

}
