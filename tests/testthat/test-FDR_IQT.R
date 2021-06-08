library(winnerscurse)

context("FDR IQT")


test_that("testing if FDR-IQT function adjusts correctly",

          {
            summary_data <- data.frame(rsid = c(1,2),beta = c(0.645,0.49), se =c(0.25,0.25))
            out <- FDR_IQT(summary_data)

            testdf <- data.frame(beta_FIQT = c(2.33*0.25,1.96*0.25)) # calculated by hand

            expect_true(identical(round(out$beta_FIQT,2),round(testdf$beta_FIQT,2))==1)

          })
