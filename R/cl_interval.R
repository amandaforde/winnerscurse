#' Confidence interval associated with discovery GWAS conditional likelihood
#' methods
#'
#' \code{cl_interval} is a function that allows the user to obtain a confidence
#' interval for the adjusted association estimates of significant SNPs, which
#' have been obtained through the implementation of
#' \code{\link{conditional_likelihood}}. This function produces one confidence
#' interval for each significant SNP, based on the approach suggested in
#' \href{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2665019/}{Ghosh \emph{et
#' al.} (2008)}. Note that in order for an appropriate confidence interval to be
#' outputted for each significant SNP, the absolute value of the largest
#' \eqn{z}-statistic in the data set must be less than 150.
#'
#' @param summary_data A data frame containing summary statistics from the
#'   discovery GWAS. It must have three columns with column names \code{rsid},
#'   \code{beta} and \code{se}, respectively, and columns \code{beta} and
#'   \code{se} must contain numerical values. Each row must correspond to a
#'   unique SNP, identified by \code{rsid}.
#' @param alpha A numerical value which specifies the desired genome-wide
#'   significance threshold. The default is given as \code{5e-8}.
#' @param conf_level A numerical value between 0 and 1 which determines the
#'   confidence interval to be computed. The default setting is \code{0.95}
#'   which results in the calculation of a 95\% confidence interval for the
#'   adjusted association estimate for each SNP.
#'
#' @return A data frame which combines the output of
#'   \code{conditional_likelihood} with two additional columns, namely
#'   \code{lower} and \code{upper}, containing the lower and upper bounds of the
#'   required confidence interval for each significant SNP, respectively.
#'   However, if no SNPs are detected as significant in the data set,
#'   \code{cl_interval} returns a warning message: \code{"WARNING: There are no
#'   significant SNPs at this threshold."}
#'
#' @references Ghosh, A., Zou, F., & Wright, F. A. (2008). Estimating odds
#'   ratios in genome scans: an approximate conditional likelihood approach.
#'   \emph{American journal of human genetics}, \strong{82(5)}, 1064\eqn{-}1074.
#'   \doi{10.1016/j.ajhg.2008.03.002}
#'
#' @seealso \code{\link{conditional_likelihood}} for details on operation of
#'   conditional likelihood methods with summary statistics from discovery GWAS.
#'
#'   \url{https://amandaforde.github.io/winnerscurse/articles/standard_errors_confidence_intervals.html}
#'    for illustration of the use of \code{cl_interval} with a toy data set and
#'   further information regarding the manner in which the confidence interval
#'   is computed.
#' @export
#'
cl_interval <- function(summary_data,alpha=5e-8, conf_level=0.95){

  stopifnot(conf_level < 1 && conf_level > 0)

  p_val <- 2*(stats::pnorm(abs(summary_data$beta/summary_data$se), lower.tail=FALSE))

  if (sum(p_val<alpha) == 0) {
    warn <- "WARNING: There are no significant SNPs at this threshold."
    return(warn)
  }

  summary_data <- conditional_likelihood(summary_data,alpha)

  z <- summary_data$beta/summary_data$se
  c <- stats::qnorm((alpha)/2, lower.tail=FALSE)

  lower <- c(rep(0),nrow(summary_data))
  upper <- c(rep(0),nrow(summary_data))

  for (i in 1:nrow(summary_data)){
    cond.like <- function (mu) { return ((stats::dnorm(z[i]-mu))/(stats::pnorm(-c+mu)+stats::pnorm(-c-mu)))}

    if ( abs(z[i]) < 37){
      norm.cond.like <- function (x) {cond.like(x)/(stats::integrate(cond.like,-37,37)$value)}

      fun_low <- function(y){return(stats::integrate(norm.cond.like,-37,y)$value - ((1-conf_level)/2))}
      lower[i] <- (stats::uniroot(fun_low,interval=c(-37,37))$root)*summary_data$se[i]

      fun_upp <- function(y){return(stats::integrate(norm.cond.like,-37,y)$value - (1-(1-conf_level)/2))}
      upper[i] <- (stats::uniroot(fun_upp,interval=c(-37,37))$root)*summary_data$se[i]

    }else{
      norm.cond.like <- function (x) {cond.like(x)/(stats::integrate(cond.like,-150,150)$value)}

      fun_low <- function(y){return(stats::integrate(norm.cond.like,-150,y)$value - ((1-conf_level)/2))}
      lower[i] <- (stats::uniroot(fun_low,interval=c(-150,150))$root)*summary_data$se[i]

      fun_upp <- function(y){return(stats::integrate(norm.cond.like,-150,y)$value - (1-(1-conf_level)/2))}
      upper[i] <- (stats::uniroot(fun_upp,interval=c(-150,150))$root)*summary_data$se[i]

    }


  }

  summary_data <- cbind(summary_data,lower,upper)
  return(summary_data)
}
