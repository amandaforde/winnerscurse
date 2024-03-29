% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/empirical_bayes.R
\name{empirical_bayes}
\alias{empirical_bayes}
\title{Empirical Bayes method for use with discovery GWAS}
\usage{
empirical_bayes(summary_data, method = "AIC")
}
\arguments{
\item{summary_data}{A data frame containing summary statistics from the
discovery GWAS. It must have three columns with column names \code{rsid},
\code{beta} and \code{se}, respectively, and columns \code{beta} and
\code{se} must contain numerical values. Each row must correspond to a
unique SNP, identified by \code{rsid}. The function requires that there must
be at least 50 SNPs as any less will result in very poor performance of the
method.}

\item{method}{A string which allows the user to choose what modelling approach
to take for the purpose of estimating the log density function. The default
setting is \code{method="AIC"}, which is the current published method. Other
options include \code{method="fix_df"}, \code{method="scam"},
\code{method="gam_nb"} and \code{method="gam_po"}. If
\code{method="fix_df"}, the degrees of freedom is set to 7. The other three
options all enforce additional constraints on the shape of the estimated log
density function.}
}
\value{
A data frame with the inputted summary data occupying the first three
 columns. The new adjusted association estimates for each SNP are returned in
 the fourth column, namely \code{beta_EB}. The SNPs are contained in this
 data frame according to their significance, with the most significant SNP,
 i.e. the SNP with the largest absolute \eqn{z}-statistic, now located in the
 first row of the data frame.
}
\description{
\code{empirical_bayes} is a function which uses summary statistics to correct
for bias induced by Winner's Curse in SNP-trait association estimates,
obtained from a discovery GWAS. The function is strongly based on the method
originally detailed in
\href{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4048064/}{Ferguson \emph{et
al.} (2013)}. However, the function also includes all potential adaptations to the empirical Bayes method
discussed in \href{https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1010546}{Forde \emph{et al.} (2023)}.
}
\references{
Ferguson, J. P., Cho, J. H., Yang, C., & Zhao, H. (2013).
 Empirical Bayes correction for the Winner's Curse in genetic association
 studies. \emph{Genetic epidemiology}, \strong{37(1)}, 60\eqn{-}68.

Forde, A., Hemani, G., & Ferguson, J. (2023).
 Review and further developments in statistical corrections for Winner’s Curse in genetic association studies. \emph{PLoS Genetics}, \strong{19(9)}, e1010546.
}
\seealso{
\url{https://amandaforde.github.io/winnerscurse/articles/winners_curse_methods.html}
for illustration of the use of \code{empirical_bayes} with a toy data set and
further information regarding the computation of the adjusted SNP-trait
association estimates.
}
