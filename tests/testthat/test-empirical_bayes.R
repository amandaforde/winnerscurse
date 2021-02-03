library(winnerscurse)

context("empirical bayes")



test_that("testing example",

          {

            summary_data <- data.frame(rsID = c(1,2),beta = c(0.49,0.645), se =c(0.25,0.25))

            out <- empirical_bayes(summary_data)

            expect_true(out == 1)

          })
