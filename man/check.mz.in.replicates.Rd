\name{check.mz.in.replicates}
\alias{check.mz.in.replicates}
\title{check.mz.in.replicates}
\description{Function to find metabolic characteristics of individuals.}
\usage{
check.mz.in.replicates(dataA, min.samps=2, min.reps=2, num_replicates=3)
}
\arguments{
  \item{dataA}{Sample intensity matrix from apLCMS or XCMS. 
  This should only include columns corresponding to samples. All other
  information such as m/z, retention time, etc. should be deleted.}
  \item{min.samps}{min.samps: minimum number of samples for which a feature signal should be detected in at least min.reps replicates~~}
  \item{min.reps}{minimum proportion of replicates in which a signal is present (eg: 0.5 or 1)}
  \item{num_replicates}{number of replicats for each sample (eg: 2)}
 }
\details{This utility allows identification of rare features that are present 
in only some biological samples, but are present in majority of the analytical
replicates of individual samples as a result of unique environmental exposure. 
The min.samps and min.reps are user defined values for defining the minimum
number of samples and minimum proportion of replicates in which a feature
should be detected.}
\value{Filtered matrix of features with sample intensities}
\author{Karan Uppal <kuppal2@emory.edu>}
\keyword{ ~replicates }
