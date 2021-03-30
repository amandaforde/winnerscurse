#' Adaptation of BR-squared method - Faye et al. (2011) - suitable for summary
#' statistics
#'
#' @param summary_data Data frame containing summary data, three columns: rsid,
#'   beta, se
#' @param seed Allows reproducibility of adjusted estimates as variances can
#'   occur due to bootstrap, default is seed=1998
#'
#' @return Data frame with summary data together with adjusted estimates
#' @export

BR_ss <- function(summary_data,seed=1998){
    summary_data <- dplyr::arrange(summary_data, dplyr::desc((summary_data$beta/summary_data$se)))

    N <- nrow(summary_data)
    set.seed(seed)
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
