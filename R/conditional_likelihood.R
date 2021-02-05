#' Conditional likelihood methods to correct for Winner’s Curse - Ghosh et al. (2008)
#'
#' @param summary_data Data frame containing summary data, three columns: rsid, beta, se
#' @param alpha significance threshold used in discovery GWAS
#'
#' @return Data frame with summary data of only significant SNPs together with three corrected estimates
#' or data frame with all summary data and corrected estimates of significant SNPs
#' @export
#'
#'
#'
conditional_likelihood <- function(summary_data, alpha){

  z <- summary_data$beta/summary_data$se
  p_val <- 2*(1-stats::pnorm(abs(z)))
  summary_data <- cbind(summary_data, z, p_val)


  if (sum(summary_data$p_val<alpha) == 0) {
    summary_data <- dplyr::arrange(summary_data,dplyr::desc(abs(summary_data$z)))
    return(summary_data[,1:3])
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
    beta.cl2[i] <- ((stats::integrate(mean.cond.like,-37,37))$value/(stats::integrate(cond.like,-37,37))$value)*(summary_data_sig$se[i])
    beta.cl3[i] <- (beta.cl1[i]+beta.cl2[i])/2
  }


  summary_data_sig <- cbind(summary_data_sig,beta.cl1,beta.cl2,beta.cl3)
  summary_data_sig <- dplyr::arrange(summary_data_sig,dplyr::desc(abs(summary_data_sig$z)))
  summary_data_sig <- cbind(summary_data_sig[,1:3],summary_data_sig[,6:8])
  return(summary_data_sig)

}
