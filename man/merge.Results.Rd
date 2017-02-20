\name{merge.Results}
\alias{merge.Results}
\title{merge.Results}
\description{ Function that merges results from different parameter settings.}
\usage{
merge.Results(dataA, dataB, feat.eval.A, feat.eval.B, max.mz.diff = 15, 
max.rt.diff = 300,merge.eval.pvalue = 0.05, alignment.tool="apLCMS", numnodes = 1, 
mult.test.cor=FALSE,mergecorthresh=0.7, missingvalue=0)
}
\arguments{
  \item{dataA}{feature alignment output matrix from apLCMS or XCMS with
 intensities at parameter settings P1}
 \item{dataB}{feature alignment output matrix from apLCMS or XCMS with
 intensities at parameter settings P2}
 \item{feat.eval.A}{feature evaluations results from evaluate.Features
 for results at parameter settings P1}
  \item{feat.eval.B}{feature evaluations results from evaluate.Features
 for results at parameter settings P2}
  \item{max.mz.diff}{+/- mz tolerance in ppm for feature matching}
  \item{max.rt.diff}{retention time tolerance for feature matching}
  \item{merge.eval.pvalue}{Threshold for defining signifcance level of the paired t-test or
  the Pearson correlation during the merge stage in xMSanalyzer. The p-value is used
  to determine whether two features with same m/z and retention time have identical
  intensity profiles.}
  \item{mergecorthresh}{Correlation threshold to be used during the merge stage in xMSanalyzer
  to determine whether two features with same m/z and retention time have identical
  intensity profiles.}
  \item{alignment.tool}{Name of the feature alignment tool eg: "apLCMS" or "XCMS"}
  \item{numnodes}{Number of computing nodes to use. Default: 1}
  \item{mult.test.cor}{Should Bonferroni multiple testing correction method 
  be applied for comparing intensities of overlapping m/z? Default: FALSE} 
  \item{missingvalue}{How are missing values represented. eg: 0 or NA}
  }
\details{ We use a four-step process to merge features from different parameter settings. 
In step one, features detected at settings P1 and P2 are combined into one list. In 
step two, features are grouped by a user-defined m/z tolerance (5 ppm is appropriate
for high resolution MS but may not be suitable for lower resolution instruments; 
for the LTQ-FT/MS, examination of m/z tolerance shows little difference between 5 and 10 ppm).
In step three, features are further sub-grouped based on a user-defined retention time tolerance.
Users are recommended to use the find.Overlapping.mzs function below to optimize the retention
time tolerance threshold. In step four, a paired t-test & a Spearman correlation test is used to 
compare the intensity levels of the metabolites only for the redundant features that have m/z 
and retention time within defined tolerance levels as described above.
Features with minimum median PID
(or median CV; for more than two technical replicates) are chosen as representatives of
each sub-group, and added to the final list. This scheme allows identification of unique
features, and selection of the most consistent feature as a representative for features
that overlap.}
\value{Returns a matrix of columns of unique m/z values, elution times, signal strengths in each spectrum
after merging results}
\author{Karan Uppal <kuppal2@emory.edu>}
\keyword{ ~merge }
