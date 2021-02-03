#' Conditional likelihood methods to correct for Winner’s Curse - Ghosh et al. (2008)
#'
#' @param summary_data Data frame containing summary data of only significant SNPs, three columns: rsid, beta, se
#' @param alpha significance threshold used in discovery GWAS
#'
#' @return Data frame with summary data together with three corrected estimates
#' @export
#'
#'
#'
conditional_likelihood <- function(summary_data, alpha)

{
  z <- summary_data$beta/summary_data$se
  c <- stats::qnorm(1-(alpha)/2)

  beta.cl1 <- c(rep(0,nrow(summary_data)))
  beta.cl2 <- c(rep(0,nrow(summary_data)))
  beta.cl3 <- c(rep(0,nrow(summary_data)))

  for (i in 1:nrow(summary_data)){
    cond.like <- function (mu) { return ((stats::dnorm(z[i]-mu))/(stats::pnorm(-c+mu)+stats::pnorm(-c-mu)))}
    mean.cond.like <- function (mu) {return(mu*cond.like(mu))}

    beta.cl1[i] <- (stats::optimize(cond.like, c(0,z[i]), maximum=TRUE)$maximum)*summary_data$se[i]
    beta.cl2[i] <- ((stats::integrate(mean.cond.like,-37,37))$value/(stats::integrate(cond.like,-37,37))$value)*(summary_data$se[i])
    beta.cl3[i] <- (beta.cl1[i]+beta.cl2[i])/2
  }

  summary_data <- cbind(summary_data,beta.cl1,beta.cl2,beta.cl3)

  return(summary_data)


}
