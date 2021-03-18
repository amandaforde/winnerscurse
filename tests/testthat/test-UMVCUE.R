library(winnerscurse)

context("UMVCUE")



test_that("testing if UMVCUE gives appropriate estimates when given two data sets",

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

            summary_disc <- data.frame(rsid=seq(1,n_snps),beta=rnorm(n=n_snps,mean=true_beta,sd=se),se=se)

            n_samples_rep <- 450000
            se_rep <- 1/sqrt(2*n_samples_rep*maf*(1-maf))
            summary_rep <- data.frame(rsid=seq(1,n_snps),beta=rnorm(n=n_snps,mean=true_beta,sd=se_rep),se=se_rep)

            out <- UMVCUE(summary_disc, summary_rep, alpha=5e-8)

            test <- sum(abs(round(out$disc_beta,6)) >= abs(round(out$beta_UMVCUE,6))) >= 0.75*length(out$disc_beta)

            expect_true(identical(test,TRUE) == 1)

          })


