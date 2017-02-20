\name{evaluate.Samples}
\alias{evaluate.Samples}
\title{evaluate.Samples}
\description{Evaluate sample consistency based on Pearson or Spearman Correlation.}
\usage{evaluate.Samples(curdata, numreplicates, alignment.tool, cormethod = "pearson",
missingvalue = 0, ignore.missing = TRUE, replace.bad.replicates = TRUE)}
\arguments{
  \item{curdata}{feature alignment output matrix from apLCMS or XCMS with intensities}
  \item{numreplicates}{number of technical replicates per sample}
  \item{alignment.tool}{name of the feature alignment tool eg: "apLCMS" or "XCMS"}
  \item{cormethod}{Pearson or Spearman correlation.}
   \item{missingvalue}{How are missing values represented? eg: 0 or NA}
  \item{ignore.missing}{Should the missing values be ignored while computing
  pearson correlation? eg: TRUE or FALSE}
  \item{replace.bad.replicates}{Should the bad replicates be replaced by the
  average of the good ones? For example, if the number of technical replicates is
  more than two, and one of the replicates is poorly correlated with the other
  two but the other replicates have correlation greater than the defined
  threshold, then the bad replicate is replaced by the average of the good ones.}
}
\details{If at least two analytical replicates are present for each biological sample,
this function calculates the mean pairwise Pearson correlation coefficient between
sample replicates using the built-in cor() function in R. Only the features with no
missing values are used to evaluate correlation. Analytical replicates refer to
multiple injections from the same biological sample; whereas, samples refer to
different biological samples.
}
\value{returns a matrix of Pearson or Spearman Correlation Coefficients within technical replicates per sample.}
\author{Karan Uppal <kuppal2@emory.edu>}
\keyword{ ~Sample quality }
