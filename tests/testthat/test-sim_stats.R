library(winnerscurse)

context("sim_stats")


test_that("testing if sim_stats function simulates summary statistics appropriately",

          {

            data <- sim_stats() ## using all default params

            expect_true(identical(is.null(data$rep),TRUE)==1)
            expect_true(identical(nrow(data$true),nrow(data$disc),10^6)==1)

            data1 <- sim_stats(rep=TRUE)

            expect_true(identical(is.null(data1$rep),FALSE)==1)
            expect_true(identical(ncol(data1$disc),ncol(data1$rep),3)==1)


          })


