\name{getCVreplicates}
\alias{getCVreplicates}
\title{getCVreplicates}
\description{Evaluate feature consistency based on coefficient of variation, where cv=sd/mean. 
Only calculate CV for samples for which a signal is detected in at least 60% of technical replicates.  
If signal is detected in 60% of technical replicates, then the missing values are replaced by mean intensity
which is calculated using non-missing replicates.}
\usage{getCVreplicates(curdata, alignment.tool, numreplicates,
min.samp.percent=0.6, impute.bool=TRUE, missingvalue)}
\arguments{
  \item{curdata}{Feature alignment output matrix from apLCMS or XCMS with sample intensities}
  \item{alignment.tool}{Name of the feature alignment tool eg: "apLCMS" or "XCMS"}
  \item{numreplicates}{Number of replicates per sample}
 \item{min.samp.percent}{If signal is detected in x proportion of technical replicates, then
 the missing values are replaced by mean intensity which is calculated using non-missing
 replicates. eg: 0.7}
  \item{impute.bool}{Should the missingvalues be replaced by mean
  of the other replicates? eg: TRUE or FALSE}
  \item{missingvalue}{How are the missing values represented? eg: 0 or NA}

}
\value{Matrix of feature consistency based on CV with columns:
mz: m/z of the feature
minCV: minimum CV between technical replicates across all samples
first_quartileCV: 25th percentile CV 
medianCV: 50th percentile CV
meanCV: average of cofficient of variations between technical replicates per sample across all samples
third_quartileCV: 75th percentile CV 
maxCV: maximum CV between technical replicates across all samples  
}
\author{Karan Uppal <kuppal2@emory.edu>}
\keyword{ ~CV }