#' Conditional likelihood methods to correct for Winner’s Curse - Ghosh et al. (2008)
#'
#' @param summary_data Data frame containing summary data, three columns: rsid, beta, se
#' @param alpha significance threshold used in discovery GWAS
#' @param sig select true to only return data frame with significant SNPs and corrected estimates
#'
#' @return Data frame with summary data of only significant SNPs together with three corrected estimates
#' @export
#'
#'
#'
conditional_likelihood <- function(summary_data, alpha, sig = FALSE)

{
  summary_data_sig <- summary_data[(2*(1-stats::pnorm(abs(summary_data$beta/summary_data$se))))<alpha,]
  summary_data_notsig <- summary_data[(2*(1-stats::pnorm(abs(summary_data$beta/summary_data$se))))>=alpha,]

  z <- summary_data_sig$beta/summary_data_sig$se
  c <- stats::qnorm(1-(alpha)/2)

  beta.cl1 <- c(rep(0,nrow(summary_data_sig)))
  beta.cl2 <- c(rep(0,nrow(summary_data_sig)))
  beta.cl3 <- c(rep(0,nrow(summary_data_sig)))

  for (i in 1:nrow(summary_data_sig)){
    cond.like <- function (mu) { return ((stats::dnorm(z[i]-mu))/(stats::pnorm(-c+mu)+stats::pnorm(-c-mu)))}
    mean.cond.like <- function (mu) {return(mu*cond.like(mu))}

    beta.cl1[i] <- (stats::optimize(cond.like, c(0,z[i]), maximum=TRUE)$maximum)*summary_data_sig$se[i]
    beta.cl2[i] <- ((stats::integrate(mean.cond.like,-37,37))$value/(stats::integrate(cond.like,-37,37))$value)*(summary_data_sig$se[i])
    beta.cl3[i] <- (beta.cl1[i]+beta.cl2[i])/2
  }

  if (sig == TRUE) {
    summary_data_sig <- cbind(summary_data_sig,beta.cl1,beta.cl2,beta.cl3)
    summary_data_sig <- dplyr::arrange(summary_data_sig,dplyr::desc(abs(summary_data_sig$beta/summary_data_sig$se)))
    return(summary_data_sig)
  }


  else {
    summary_data_sig <- cbind(summary_data_sig,beta.cl1,beta.cl2,beta.cl3)
    beta.cl_notsig <- data.frame(beta.cl1 = summary_data_notsig$beta, beta.cl2 = summary_data_notsig$beta, beta.cl3 = summary_data_notsig$beta)
    summary_data_notsig <- cbind(summary_data_notsig,beta.cl_notsig)
    summary_data <- rbind(summary_data_sig,summary_data_notsig)
    summary_data <- dplyr::arrange(summary_data,dplyr::desc(abs(summary_data$beta/summary_data$se)))
    return(summary_data)
  }

}
