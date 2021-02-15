library(winnerscurse)

context("BR ss")



test_that("testing if BR_ss is giving appropriate output",

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

            out <- BR_ss(summary_stats)
            out_sig <- out[2*(1-pnorm(abs(out$beta/out$se))) < 5e-8,]

            test <- sum(abs(round(out_sig$beta,6)) >= abs(round(out_sig$beta_BR_ss,6))) == length(out_sig$beta)

            expect_true(identical(test,TRUE) == 1)

          })
