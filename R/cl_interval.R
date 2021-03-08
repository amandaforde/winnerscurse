#' Confidence interval for conditional likelihood methods - Ghosh et al. (2008)
#'
#' @param summary_data Data frame containing summary data, three columns: rsid, beta, se
#' @param alpha significance threshold used in discovery GWAS
#' @param conf specifies confidence level, takes a value between 0 and 1 - default is 0.95 corresponding to 95 per cent confidence interval
#'
#' @return Data frame which combines output from application of conditional likelihood methods, i.e. summary data of only significant SNPs together with three corrected estimates, with an upper and a lower bound on the required confidence interval for each SNP
#' @export
#'
cl_interval <- function(summary_data,alpha, conf=0.95){
  summary_data <- conditional_likelihood(summary_data,alpha)

  z <- summary_data$beta/summary_data$se
  c <- stats::qnorm(1-(alpha)/2)

  lower <- c(rep(0),nrow(summary_data))
  upper <- c(rep(0),nrow(summary_data))

  for (i in 1:nrow(summary_data)){
    cond.like <- function (mu) { return ((stats::dnorm(z[i]-mu))/(stats::pnorm(-c+mu)+stats::pnorm(-c-mu)))}
    norm.cond.like <- function (x) {cond.like(x)/(stats::integrate(cond.like,-37,37)$value)}

    fun_low <- function(y){return(stats::integrate(norm.cond.like,-37,y)$value - ((1-conf)/2))}
    lower[i] <- (stats::uniroot(fun_low,interval=c(-37,37))$root)*summary_data$se[i]

    fun_upp <- function(y){return(stats::integrate(norm.cond.like,-37,y)$value - (1-(1-conf)/2))}
    upper[i] <- (stats::uniroot(fun_upp,interval=c(-37,37))$root)*summary_data$se[i]

  }

  summary_data <- cbind(summary_data,lower,upper)
  return(summary_data)
}
