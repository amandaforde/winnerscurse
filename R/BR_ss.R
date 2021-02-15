#' Adaptation of BR-squared method - Faye et al. (2011) - suitable for summary statistics
#'
#' @param summary_data Data frame containing summary data, three columns: rsid, beta, se
#'
#' @return Data frame with summary data together with adjusted estimates
#' @export

BR_ss <- function(summary_data){

  summary_data <- dplyr::arrange(summary_data, dplyr::desc(abs(summary_data$beta/summary_data$se)))

  N <- nrow(summary_data)
  beta_boot <- matrix(stats::rnorm(100*N, mean = rep(summary_data$beta,100), sd = rep(summary_data$se,100)), nrow=N, ncol=100, byrow=FALSE)

  beta_mat <- matrix(rep(summary_data$beta,100), nrow=N, ncol=100, byrow=FALSE)
  se_mat <- matrix(rep(summary_data$se,100), nrow=N, ncol=100, byrow=FALSE)

  beta_oob <- beta_mat

  ordering <- apply(beta_boot/se_mat, 2,order,decreasing=TRUE)
  bias_correct <- matrix(nrow=N, ncol=100)

  for (i in 1:100){
    bias_correct[,i] <- (beta_boot[ordering[,i],i] - beta_oob[ordering[,i],i])/summary_data$se[ordering[,i]]
  }
  bias_correct <- apply(bias_correct,1,mean)

  beta_BR_ss <- summary_data$beta - summary_data$se*bias_correct[rank(-1*summary_data$beta/summary_data$se)]
  beta_BR_ss[sign(beta_BR_ss) != sign(summary_data$beta)] <- 0

  summary_data <- cbind(summary_data, beta_BR_ss)

  return(summary_data)
}
