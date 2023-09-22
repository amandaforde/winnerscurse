library(winnerscurse)

context("conditional likelihood confidence interval")



test_that("testing if conditional likelihood confidence intervals outputted correctly",

          {
            summary_data <- data.frame(rsid = c(1,2),beta = c(1.5,1.3), se =c(0.25,0.25))
            alpha <- 2*(1-pnorm(5))

            out <- cl_interval(summary_data,alpha)
            testdf <- data.frame(beta.cl1 = c(5.48,0.66),beta.cl2 = c(4.94,2.53),beta.cl3 = c(5.21,1.6))
            # above results were quoted in Ghosh et al. (2008)
            expect_true(identical(round(out[,4:6]/out[,3],2),testdf)==1)

            test <- sum(out$lower < out$upper) == length(out$rsid)
            expect_true(identical(test,TRUE)==1)

            out <- cl_interval(summary_data,alpha=1e-10)
            test1 <- out == "WARNING: There are no significant SNPs at this threshold."
            expect_true(identical(test1,TRUE)==1)

          })

