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

            test <- sum(abs(round(out_sig$beta,6)) >= abs(round(out_sig$beta_BR_ss,6))) >= length(out_sig$beta)*0.7
            expect_true(identical(test,TRUE) == 1)


            maf <- c(0.22,0.42,0.42,0.08,0.39,0.29,0.47,0.49,0.14,0.13)
            se <- 1/sqrt(2*30000*maf*(1-maf))
            beta <- c(0.08,0.59,0.1,0.1,0.3,-0.1,-0.95,0.04,-0.29,0.05)
            summary_stats <- data.frame(rsid=seq(1,10),beta=beta,se=se)

            out1 <- BR_ss(summary_stats,seed_opt=TRUE)
            out2 <- BR_ss(summary_stats,seed_opt=TRUE)

            test <- as.numeric(sum(out1 == out2))
            expect_true(identical(test,40) == 1)

          })
