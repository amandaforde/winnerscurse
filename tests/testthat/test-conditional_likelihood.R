library(winnerscurse)

context("conditional likelihood")



test_that("testing if conditional likelihood functions adjust correctly",

          {
            summary_data <- data.frame(rsID = c(1,2),beta = c(1.3,1.5), se =c(0.25,0.25))
            alpha <- 2*(1-pnorm(5))

            out <- conditional_likelihood(summary_data,alpha)
            testdf <- data.frame(beta.cl1 = c(0.66,5.48),beta.cl2 = c(2.53,4.94),beta.cl3 = c(1.6,5.21))
            # above results were quoted in Ghosh et al. (2008)

            expect_true(identical(round(out[,4:6]/out[,3],2),testdf)==1)

          })
