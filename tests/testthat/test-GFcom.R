library(winnerscurse)

context("GFcom")



test_that("testing if GFcom gives appropriate estimates when given two data sets",

          {
            set.seed(1998)
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

            n_samples_rep <- 45000
            se_rep <- 1/sqrt(2*n_samples_rep*maf*(1-maf))
            summary_rep <- data.frame(rsid=seq(1,n_snps),beta=rnorm(n=n_snps,mean=true_beta,sd=se_rep),se=se_rep)

            out <- GFcom(summary_disc, summary_rep, alpha=5e-8)

            test <- sum(round(max(abs(out$beta_disc),abs(out$beta_rep)),6) >= abs(round(out$beta_GFcom,6))) == length(out$beta_disc)

            expect_true(identical(test,TRUE) == 1)


            summary_disc_1 <- summary_disc[abs(summary_disc$beta/summary_disc$se) < stats::qnorm((5e-8)/2, lower.tail=FALSE),]
            summary_rep_1 <- summary_disc[abs(summary_disc$beta/summary_disc$se) < stats::qnorm((5e-8)/2, lower.tail=FALSE),]

            out <- GFcom(summary_disc_1, summary_rep_1, alpha=5e-8)
            test2 <- sum(colnames(out) == c("rsid","beta_disc","se_disc","beta_rep","se_rep")) == 5
            expect_true(identical(test2,TRUE) == 1)


            sig1 <- which(abs(summary_disc$beta/summary_disc$se) > stats::qnorm((5e-8)/2, lower.tail=FALSE))[1]
            summary_disc_2 <- rbind(summary_disc_1,summary_disc[sig1,])
            summary_rep_2 <-  rbind(summary_rep_1,summary_rep[sig1,])

            out <- GFcom(summary_disc_2, summary_rep_2, alpha=5e-8)
            test3 <- nrow(out) == 1 && ncol(out) == 6
            expect_true(identical(test3,TRUE) == 1)


          })

