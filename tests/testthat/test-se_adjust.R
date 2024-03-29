library(winnerscurse)

context("se adjust")



test_that("testing if se_adjust gives bootstrap standard errors of adjusted estimates",

          {
            n_snps <- 10^6
            effect_snps <- 10000
            n_samples <- 30000

            maf <- runif(n_snps,0.01,0.5)
            se <- 1/sqrt(2*n_samples*maf*(1-maf))

            true_beta <- rnorm(effect_snps,0,1)
            h2 <- 0.7 # variance explained by effect SNPs
            var_y <- sum(2*maf[1:effect_snps]*(1-maf[1:effect_snps])*true_beta^2)/h2

            true_beta <- true_beta/sqrt(var_y) # scaling to represent a phenotype with variance 1
            true_beta <- c(true_beta, rep(0,n_snps-effect_snps))

            summary_stats <- data.frame(rsid=seq(1,n_snps),beta=rnorm(n=n_snps,mean=true_beta,sd=se),se=se)

            out <- se_adjust(summary_stats, method="empirical_bayes", n_boot=10)

            test <- sum(abs(round(out$beta,6)) >= abs(round(out$beta_EB,6))) >= 0.9*length(out$beta)
            test2 <- sum(out$adj_se <= out$se) >= 0.9*length(out$se) ## behaviour observed previously!

            expect_true(identical(test,TRUE) == 1)
            expect_true(identical(test2,TRUE) == 1)


            out <- se_adjust(summary_stats, method="FDR_IQT", n_boot=10)

            test3 <- sum(abs(round(out$beta,6)) >= abs(round(out$beta_FIQT,6))) >= 0.9*length(out$beta)
            test4 <- sum(out$adj_se <= out$se) >= 0.5*length(out$se)

            expect_true(identical(test3,TRUE) == 1)
            expect_true(identical(test4,TRUE) == 1)

            out <- se_adjust(summary_stats, method="BR_ss", n_boot=10)

            test5 <- sum(abs(round(out$beta,6)) >= abs(round(out$beta_BR_ss,6))) >= 0.9*length(out$beta)
            test6 <- sum(out$adj_se <= out$se) >= 0.5*length(out$se)

            expect_true(identical(test5,TRUE) == 1)
            expect_true(identical(test6,TRUE) == 1)

          })

