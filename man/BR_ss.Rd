% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BR_ss.R
\name{BR_ss}
\alias{BR_ss}
\title{Bootstrap method for use with discovery GWAS}
\usage{
BR_ss(summary_data, seed_opt = FALSE, seed = 1998)
}
\arguments{
\item{summary_data}{A data frame containing summary statistics from the
discovery GWAS. It must have three columns with column names \code{rsid},
\code{beta} and \code{se}, respectively, and all columns must contain
numerical values. Each row must correspond to a unique SNP, identified by
the numerical value \code{rsid}. The function requires that there must be at
least 5 SNPs as any less will result in issues upon usage of the
smoothing spline.}

\item{seed_opt}{A logical value which allows the user to choose if they wish
to set a seed, in order to ensure reproducibility of adjusted estimates.
Small differences can occur between iterations of the function with the same
data set due to the use of parametric bootstrapping. The default setting is
\code{seed_opt=FALSE}.}

\item{seed}{A numerical value which specifies the seed used if
\code{seed_opt=TRUE}. The default setting is the arbitrary value of
\code{1998}.}
}
\value{
A data frame with the inputted summary data occupying the first three
 columns. The new adjusted association estimates for each SNP are returned in
 the fourth column, namely \code{beta_BR_ss}. The SNPs are contained in this
 data frame according to their significance, with the most significant SNP,
 i.e. the SNP with the largest absolute \eqn{z}-statistic, now located in the
 first row of the data frame.
}
\description{
\code{BR_ss} is a function which aims to use summary statistics to alleviate
Winner's Curse bias in SNP-trait association estimates, obtained from a
discovery GWAS. The function implements an adaptation of a bootstrap
resampling method known as BR-squared, detailed in
\href{https://onlinelibrary.wiley.com/doi/abs/10.1002/sim.4228}{Faye \emph{et
al.} (2011)}.
}
\references{
Faye, L. L., Sun, L., Dimitromanolakis, A., & Bull, S. B. (2011).
 A flexible genome-wide bootstrap method that accounts for ranking and
 threshold-selection bias in GWAS interpretation and replication study
 design. \emph{Statistics in Medicine}, \strong{30(15)}, 1898\eqn{-}1912.
 \url{https://doi.org/10.1002/sim.4228}
}
\seealso{
\url{https://amandaforde.github.io/winnerscurse/articles/winners_curse_methods.html}
for illustration of the use of \code{BR_ss} with a toy data set and further
information regarding the computation of the adjusted SNP-trait association
estimates.
}
