#'GFcom method for use with discovery and replication GWASs
#'
#'\code{GFcom} is a function which attempts to produce less biased SNP-trait
#'association estimates for SNPs deemed significant in the discovery GWAS, using
#'summary statistics from both discovery and replication GWASs. The function
#'implements the method described in
#'\href{https://journals.sagepub.com/doi/epub/10.1177/0962280217741771}{Hu et
#'al. (2017)}.
#'
#'@param summary_disc A data frame containing summary statistics from the
#'  \emph{discovery} GWAS. It must have three columns with column names
#'  \code{rsid}, \code{beta} and \code{se}, respectively, and columns
#'  \code{beta} and \code{se} must contain numerical values. Each row must
#'  correspond to a unique SNP, identified by \code{rsid}.
#'@param summary_rep A data frame containing summary statistics from the
#'  \emph{replication} GWAS. It must have three columns with column names
#'  \code{rsid}, \code{beta} and \code{se}, respectively, and all columns must
#'  contain numerical values. Each row must correspond to a unique SNP,
#'  identified by the numerical value \code{rsid}. SNPs must be ordered in the
#'  exact same manner as those in \code{summary_disc}, i.e.
#'  \code{summary_rep$rsid} must be equivalent to \code{summary_disc$rsid}.
#'@param alpha A numerical value which specifies the desired genome-wide
#'  significance threshold for the discovery GWAS. The default is given as
#'  \code{5e-8}.
#'
#'@return A data frame with summary statistics and adjusted association
#'  estimates of only those SNPs which have been deemed significant in the
#'  discovery GWAS according to the specified threshold, \code{alpha}, i.e. SNPs
#'  with \eqn{p}-values less than \code{alpha}. The inputted summary data
#'  occupies the first five columns, in which the columns \code{beta_disc} and
#'  \code{se_disc} contain the statistics from the discovery GWAS and columns
#'  \code{beta_rep} and \code{se_rep} hold the replication GWAS statistics. The new adjusted
#'  association estimate for each SNP, as defined in the aforementioned paper,
#'  are contained in the next column, namely \code{beta_GFcom}. The SNPs are contained in this data frame according to
#'  their significance, with the most significant SNP, i.e. the SNP with the
#'  largest absolute \eqn{z}-statistic, now located in the first row of the data
#'  frame. If no SNPs are detected as significant in the discovery GWAS,
#'  \code{GFcom} merely returns a data frame which combines the two inputted
#'  data sets.
#'
#'@references Hu, J., Zhang, W., Li, X., Pan, D., & Li, Q. (2019). Efficient
#'  estimation of disease odds ratios for follow-up genetic association studies.
#'  \emph{Statistical Methods in Medical Research.}, \strong{28(7)}, 1927\eqn{-}1941.
#'  \doi{10.1177/0962280217741771}
#'
#'@export
#'
GFcom <- function(summary_disc,summary_rep,alpha=5e-8){
  ### only selection on snps from discovery data set, just require definition for c or equivalently, alpha

  stopifnot(all(c("rsid", "beta","se") %in% names(summary_disc)))
  stopifnot(!all(is.na(summary_disc$rsid)) && !all(is.na(summary_disc$beta)) && !all(is.na(summary_disc$se)))
  stopifnot(is.numeric(summary_disc$beta) && is.numeric(summary_disc$se))
  stopifnot(!any(duplicated(summary_disc$rsid)))

  stopifnot(all(c("rsid", "beta","se") %in% names(summary_rep)))
  stopifnot(!all(is.na(summary_rep$rsid)) && !all(is.na(summary_rep$beta)) && !all(is.na(summary_rep$se)))
  stopifnot(is.numeric(summary_rep$beta)  && is.numeric(summary_rep$se))
  stopifnot(!any(duplicated(summary_rep$rsid)))

  stopifnot(summary_disc$rsid == summary_rep$rsid)

  se_com <- sqrt((((summary_disc$se)^2)*((summary_rep$se)^2))/(((summary_disc$se)^2) + ((summary_rep$se)^2)))
  beta_com <- ((((summary_rep$se)^2)*(summary_disc$beta))+(((summary_disc$se)^2)*(summary_rep$beta)))/(((summary_disc$se)^2) + ((summary_rep$se)^2))
  summary_com <- data.frame(rsid=summary_disc$rsid,beta=beta_com,se=se_com)
  c_1 <- stats::qnorm((alpha)/2, lower.tail=FALSE)

  if(sum(abs(summary_disc$beta/summary_disc$se) > c_1) == 0){
    summary_data <- cbind(summary_disc[1:3],summary_rep[2:3])
    names(summary_data)[2] <- "beta_disc"
    names(summary_data)[3] <- "se_disc"
    names(summary_data)[4] <- "beta_rep"
    names(summary_data)[5] <- "se_rep"
    return(summary_data)
  }

  summary_disc_sig <- summary_disc[abs(summary_disc$beta/summary_disc$se) > c_1,]
  summary_rep_sig <- summary_rep[abs(summary_disc$beta/summary_disc$se) > c_1,]
  summary_com_sig <- summary_com[abs(summary_disc$beta/summary_disc$se) > c_1,]

  MLE.beta.com <- function(beta.com,se1,se2,C,se.com){
    L = function(beta,x){
      g_x = function(x){
        pnorm((x-C*se1)/(se.com*se1/se2))+pnorm((-x-C*se1)/(se.com*se1/se2));
      }
      numerator = dnorm((x-beta)/se.com,log = TRUE)+log(g_x(x));
      denominator.func = function(t){
        (dnorm((t-beta)/se.com)+dnorm((t+beta)/se.com))*g_x(t);
      }
      all = numerator- log(integrate(denominator.func,0,Inf)$value);
      -all;
    }
    optims = optim(par = beta.com,fn=L,x=beta.com,method = 'L-BFGS-B',lower = -1,upper=1,control = list(maxit = 1e7),hessian=TRUE);
    beta.MLE = optims$par;
    convergence = optims$convergence;
    hessian = optims$hessian[1,1];
    while(convergence!=0 | hessian<0){
      s = runif(1,-1,1);
      optims = optim(par = s,fn=L,x=beta.com,method = 'L-BFGS-B',lower = -1,upper=1,control = list(maxit = 1e7),hessian =TRUE);
      if(optims$convergence ==0 & optims$hessian[1,1]>0){
        beta.MLE = optims$par;
        hessian = optims$hessian[1,1];
      }
      convergence = optims$convergence;
    }
    se.MLE = sqrt(1/hessian);
    c(beta.MLE,se.MLE);
  }

  beta_MLE <- c(rep(0,length(summary_com_sig$beta)))
  se_MLE <- c(rep(0,length(summary_com_sig$beta)))
  for (b in 1:length(beta_MLE)){
    MLE <- MLE.beta.com(summary_com_sig$beta[b],summary_disc_sig$se[b],summary_rep_sig$se[b],c_1,summary_com_sig$se[b])
    beta_MLE[b] <- MLE[1]
    se_MLE[b] <- MLE[2]
  }

  MLE_mse <- (se_MLE)^2 + (beta_MLE - summary_rep_sig$beta)^2
  K1 <- (summary_rep_sig$se)^2/(MLE_mse + (summary_rep_sig$se)^2)
  K2 <- MLE_mse/(MLE_mse + (summary_rep_sig$se)^2)
  beta_GFcom <- K2*summary_rep_sig$beta + K1*beta_MLE

  betas <- data.frame(rsid = summary_disc_sig$rsid, beta_disc = summary_disc_sig$beta, se_disc  = summary_disc_sig$se, beta_rep = summary_rep_sig$beta, se_rep = summary_rep_sig$se, beta_GFcom=beta_GFcom)

  ## reorder based on significance in first set up!
  out <- dplyr::arrange(betas, dplyr::desc(abs(betas$beta_disc/betas$se_disc)))

  return(out)
}
