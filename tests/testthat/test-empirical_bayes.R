library(winnerscurse)

context("empirical bayes")



test_that("testing if empirical bayes is giving appropriate output",

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

            out <- empirical_bayes(summary_stats)
            test1 <- sum(abs(round(out$beta,6)) >= abs(round(out$beta_EB,6))) >= 0.9*length(out$beta)
            expect_true(identical(test1,TRUE) == 1)

            out <- empirical_bayes(summary_stats, method="fix_df")
            test2 <- sum(abs(round(out$beta,6)) >= abs(round(out$beta_EB,6))) >= 0.9*length(out$beta)
            expect_true(identical(test2,TRUE) == 1)

            out <- empirical_bayes(summary_stats, method="scam")
            test3 <- sum(abs(round(out$beta,6)) >= abs(round(out$beta_EB,6))) >= 0.9*length(out$beta)
            expect_true(identical(test3,TRUE) == 1)

            out <- empirical_bayes(summary_stats, method="gam_po")
            test4 <- sum(abs(round(out$beta,6)) >= abs(round(out$beta_EB,6))) >= 0.9*length(out$beta)
            expect_true(identical(test4,TRUE) == 1)

            out <- empirical_bayes(summary_stats, method="gam_nb")
            test5 <- sum(abs(round(out$beta,6)) >= abs(round(out$beta_EB,6))) >= 0.9*length(out$beta)
            expect_true(identical(test5,TRUE) == 1)

          })
