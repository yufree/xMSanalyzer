\name{getPID}
\alias{getPID}
\title{getPID}
\description{ Evaluate feature consistency based on PID. PID is defined as percent ratio of absolute intensity
difference to mean intensity deviation to mean intensity. Only samples with no missing values within  
technical replicates are used to evaluate PID.}
\usage{getPID(curdata, alignment.tool,missingvalue)}
\arguments{
 \item{curdata}{Feature alignment output matrix from apLCMS or XCMS with sample intensities}
  \item{alignment.tool}{Name of the feature alignment tool eg: "apLCMS" or "XCMS"}
 \item{missingvalue}{How are the missing values represented? eg: 0 or NA}

}
\value{
Matrix of feature consistency based on PID with columns:
mz: m/z of the feature
min: minimum PID between technical replicates across all samples
first_quartile: 25th percentile PID
median: 50th percentile PID
mean: average of cofficient of variations between technical replicates per sample across all samples
third_quartile: 75th percentile PID
max: maximum PID between technical replicates across all samples
}
\author{Karan Uppal <kuppal2@emory.edu>}
\keyword{ ~PID }
\keyword{ ~Feature quality }
