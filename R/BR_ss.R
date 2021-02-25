#' Adaptation of BR-squared method - Faye *et al.* (2011) - suitable for summary statistics
#'
#' @param summary_data Data frame containing summary data, three columns: rsid, beta, se
#'
#' @return Data frame with summary data together with adjusted estimates
#' @export

BR_ss <- function(summary_data){
    summary_data <- dplyr::arrange(summary_data, dplyr::desc((summary_data$beta/summary_data$se)))

    N <- nrow(summary_data)
    beta_boot <- matrix(stats::rnorm(1*N, mean = rep(summary_data$beta,1), sd = rep(summary_data$se,1)), nrow=N, ncol=1, byrow=FALSE)

    beta_mat <- matrix(rep(summary_data$beta,1), nrow=N, ncol=1, byrow=FALSE)
    se_mat <- matrix(rep(summary_data$se,1), nrow=N, ncol=1, byrow=FALSE)

    beta_oob <- beta_mat

    ordering <- apply(beta_boot/se_mat, 2,order,decreasing=TRUE)
    bias_correct <- matrix(nrow=N, ncol=1)


    bias_correct[,1] <- (beta_boot[ordering[,1],1] - beta_oob[ordering[,1],1])/summary_data$se[ordering[,1]]


    z <- summary_data$beta/summary_data$se

    index_1 <- 0.01*N
    index_2 <- N-index_1+1
    index_3 <- index_1+1
    index_4 <- N-index_1

    data_up <- data.frame(x = z[1:index_1], y = bias_correct[1:index_1])
    loess_up <- stats::predict(stats::loess(y ~ x, data = data_up, span=0.1))

    data_down <- data.frame(x = z[index_2:N], y = bias_correct[index_2:N])
    loess_down <- stats::predict(stats::loess(y ~ x, data = data_down, span=0.1))

    data <- data.frame(x = z[index_3:index_4], y = bias_correct[index_3:index_4])
    linear <- stats::predict(stats::lm(y~x,data=data))

    bias_correct <- c(loess_up,linear,loess_down)


    beta_BR_ss <- summary_data$beta - summary_data$se*bias_correct[rank(-1*summary_data$beta/summary_data$se)]
    beta_BR_ss[sign(beta_BR_ss) != sign(summary_data$beta)] <- 0

    summary_data <- cbind(summary_data, beta_BR_ss)

    summary_data <- dplyr::arrange(summary_data,dplyr::desc(abs(summary_data$beta/summary_data$se)))

    return(summary_data)
  }
