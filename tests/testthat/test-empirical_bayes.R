library(winnerscurse)

context("empirical bayes")



test_that("testing example",

          {
            N <- 10^5
            N_poly_bg <- 1000
            beta_true <- c(rnorm(N-N_poly_bg,sd=.1),rnorm(N_poly_bg,sd=1))
            MAF <- runif(N,min=0.02,max=0.5)
            SE <- 1/sqrt(2*MAF*(1-MAF))
            beta_hat <- rnorm(length(beta_true),mean=beta_true,sd=SE)

            summary_data <- data.frame(rsid = seq(1,N), beta = beta_hat, se = SE)

            out <- empirical_bayes(summary_data)

            test <- sum(abs(round(out$beta,6)) >= abs(round(out$beta_EB,6))) == length(out$beta)

            expect_true(identical(test,TRUE) == 1)

          })
