#' Empirical Bayes correction method for Winner's Curse - simplest implementation - Ferguson et al. (2013)
#'
#' @param summary_data Data frame containing summary data, three columns: rsid, beta, se
#'
#' @return Data frame with summary data together with adjusted estimates but (1) for now
#' @export
#'
#'
empirical_bayes <- function(summary_data)

{
  z <- summary_data$beta/summary_data$se

  bins <- seq(min(z),max(z),length.out=120)
  mids <- (bins[-length(bins)]+bins[-1])/2
  counts <- graphics::hist(z,breaks=bins,plot=F)$counts

  f <- stats::glm(counts ~ splines::ns(mids,df=7),stats::poisson,weight=rep(10^-50,length(counts)))$fit
  f[f==0] <- min(f[f>0])

  log_f <- as.vector(log(f))
  diff <- diff(log_f)/diff(mids)
  mids2 <- (mids[-length(mids)]+mids[-1])/2

  diff_interpol <- stats::approx(mids2,diff,mids,rule=2,ties=mean)$y

  mids_est <- c(rep(0,length(mids)))
  mids_est[mids>0] <- pmax(0, mids[mids>0] + diff_interpol[mids>0])
  mids_est[mids<0] <- pmin(0, mids[mids<0] + diff_interpol[mids<0])

  z_hat <- stats::approx(mids,mids_est,z,rule=2,ties=mean)$y
  z_hat <- sign(z)*pmin(abs(z),abs(z_hat))

  beta_EB <- z_hat*summary_data$se
  summary_data <- cbind(summary_data,beta_EB)

  return(1)

}

