#'Empirical Bayes method for use with discovery GWAS
#'
#'\code{empirical_bayes} is a function which uses summary statistics to correct
#'for bias induced by Winner's Curse in SNP-trait association estimates,
#'obtained from a discovery GWAS. The function is strongly based on the method
#'originally detailed in
#'\href{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4048064/}{Ferguson \emph{et
#'al.} (2013)}. However, the function also includes all potential adaptations to the empirical Bayes method
#'discussed in \href{https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1010546}{Forde \emph{et al.} (2023)}.
#'@param summary_data A data frame containing summary statistics from the
#'  discovery GWAS. It must have three columns with column names \code{rsid},
#'  \code{beta} and \code{se}, respectively, and columns \code{beta} and
#'  \code{se} must contain numerical values. Each row must correspond to a
#'  unique SNP, identified by \code{rsid}. The function requires that there must
#'  be at least 50 SNPs as any less will result in very poor performance of the
#'  method.
#'@param method A string which allows the user to choose what modelling approach
#'  to take for the purpose of estimating the log density function. The default
#'  setting is \code{method="AIC"}, which is the current published method. Other
#'  options include \code{method="fix_df"}, \code{method="scam"},
#'  \code{method="gam_nb"} and \code{method="gam_po"}. If
#'  \code{method="fix_df"}, the degrees of freedom is set to 7. The other three
#'  options all enforce additional constraints on the shape of the estimated log
#'  density function.
#'
#'@return A data frame with the inputted summary data occupying the first three
#'  columns. The new adjusted association estimates for each SNP are returned in
#'  the fourth column, namely \code{beta_EB}. The SNPs are contained in this
#'  data frame according to their significance, with the most significant SNP,
#'  i.e. the SNP with the largest absolute \eqn{z}-statistic, now located in the
#'  first row of the data frame.
#'
#'
#'@references Ferguson, J. P., Cho, J. H., Yang, C., & Zhao, H. (2013).
#'  Empirical Bayes correction for the Winner's Curse in genetic association
#'  studies. \emph{Genetic epidemiology}, \strong{37(1)}, 60\eqn{-}68.
#'@references Forde, A., Hemani, G., & Ferguson, J. (2023).
#'  Review and further developments in statistical corrections for Winner’s Curse in genetic association studies. \emph{PLoS Genetics}, \strong{19(9)}, e1010546.

#'
#'@seealso
#'\url{https://amandaforde.github.io/winnerscurse/articles/winners_curse_methods.html}
#'for illustration of the use of \code{empirical_bayes} with a toy data set and
#'further information regarding the computation of the adjusted SNP-trait
#'association estimates.
#'@export
#'
#'
empirical_bayes <- function(summary_data, method="AIC"){

  stopifnot(all(c("rsid", "beta","se") %in% names(summary_data)))
  stopifnot(!all(is.na(summary_data$rsid)) && !all(is.na(summary_data$beta)) && !all(is.na(summary_data$se)))
  stopifnot(is.numeric(summary_data$beta) && is.numeric(summary_data$se))
  stopifnot(nrow(summary_data) > 50)
  stopifnot(!any(duplicated(summary_data$rsid)))

  z <- summary_data$beta/summary_data$se

  bins <- seq(min(z),max(z),length.out=120)
  mids <- (bins[-length(bins)]+bins[-1])/2
  counts <- graphics::hist(z,breaks=bins,plot=F)$counts
  data <- data.frame(counts,mids)

  most_extreme <- 10
  boundary_lower <- sort(z)[most_extreme]
  boundary_upper <- sort(z,decreasing=TRUE)[most_extreme]

  if(method=="AIC"){
    df <- 7
    AIC_vector <- c(rep(0,28))
    for (best_df in 3:30){
      model <- stats::glm(counts ~ splines::ns(mids,knots = (seq(from=boundary_lower,to=boundary_upper,length=best_df+1)[2:best_df]), Boundary.knots=c(boundary_lower,boundary_upper)),stats::poisson,weights=rep(10^-50,length(counts)))
      minus2loglike <- 10^50*(model$deviance)
      AIC_vector[best_df-2] <- minus2loglike + 2*(best_df-2)
    }
    df <- 2 + which.min(AIC_vector)
    f <- stats::glm(counts ~ splines::ns(mids,knots = (seq(from=boundary_lower,to=boundary_upper,length=df+1)[2:df]), Boundary.knots=c(boundary_lower,boundary_upper)),stats::poisson,weight=rep(10^-50,length(counts)))$fit
  }

  if(method=="fix_df"){
    df <- 7
    f <- stats::glm(counts ~ splines::ns(mids,knots = (seq(from=boundary_lower,to=boundary_upper,length=df+1)[2:df]), Boundary.knots=c(boundary_lower,boundary_upper)),stats::poisson,weight=rep(10^-50,length(counts)))$fit
  }


  if(method=="gam_nb"){f <- mgcv::gam(counts ~ s(mids), family=mgcv::nb(), data=data)$fit}

  if(method=="gam_po"){f <- mgcv::gam(counts ~ s(mids), family=stats::poisson(), data=data)$fit}

  if(method=="scam"){
    f1 <- scam::scam(counts[mids >= 0] ~ s(mids[mids >= 0], bs="mpd"), family =stats::poisson(link="log"), data = data)$fit
    f2 <- scam::scam(counts[mids < 0] ~ s(mids[mids < 0], bs="mpi"), family = stats::poisson(link="log"), data = data)$fit
    f <- c(f2,f1)
  }


  f[f==0] <- min(f[f>0])
  log_f <- as.vector(log(f))
  diff <- diff(log_f)/diff(mids)
  mids2 <- (mids[-length(mids)]+mids[-1])/2

  diff_interpol <- stats::approx(mids2,diff,mids,rule=2,ties=mean)$y
  mids_est <- c(rep(0,length(mids)))
  mids_est[mids>0] <- pmax(0, mids[mids>0] + diff_interpol[mids>0])
  mids_est[mids<0] <- pmin(0, mids[mids<0] + diff_interpol[mids<0])

  z_hat <- stats::approx(mids,mids_est,z,rule=1,ties=mean)$y
  z_hat[is.na(z_hat) & z > 0] <-  z[is.na(z_hat) & z > 0] +   diff_interpol[length(diff_interpol)]
  z_hat[is.na(z_hat) & z < 0] <-  z[is.na(z_hat) & z < 0] + diff_interpol[1]
  z_hat <- sign(z)*pmin(abs(z),abs(z_hat))
  beta_EB <- z_hat*summary_data$se
  summary_data <- cbind(summary_data,beta_EB)
  summary_data <- dplyr::arrange(summary_data,dplyr::desc(abs(summary_data$beta/summary_data$se)))

  return(summary_data)

}



