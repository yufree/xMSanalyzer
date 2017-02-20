\name{evaluate.Features}
\alias{evaluate.Features}
\title{evaluate.Features}
\description{Evaluates feature consistency within analytical/technical replicates of samples based on PID or CV.
Reports quantile distribution of PID or CV within technical replicates across all samples and computes a 
quantitative reproducibility score, QRscore, which is defined as the ratio of percentage of biological samples 
for which more than 50 percent of technical replicates have a signal and median PID or CV }
\usage{
evaluate.Features(curdata, numreplicates, min.samp.percent = 0.6,alignment.tool,
impute.bool,missingvalue=0)
}
\arguments{
  \item{curdata}{Feature table from apLCMS or XCMS with sample intensities.}
  \item{numreplicates}{Number of replicates per sample. eg: 2}
  \item{min.samp.percent}{If signal is detected in x proportion of technical replicates, then the missing
  values are replaced by mean intensity which is calculated using non-missing replicates. 
  Defaul: 0.6}
  \item{alignment.tool}{Name of the feature alignment tool eg: "apLCMS" or "XCMS".}
  \item{impute.bool}{Logical (TRUE or FALSE). Used for calculating coefficient of 
	  variation (CV). Missing values are replaced by the mean of the
	  non-missing intensities if at least 60% of analytical 
	  replicates have values. For example, if two out of three 
	 replicates have non-missing values (intensity >0) and the
	third replicate has a missing value (intensity=0), then the
	   intensity of the third replicate is replaced by the mean of the 
	intensities of the other two to calculate CV.}
   \item{missingvalue}{How are the missing values represented? eg: 0 or NA}
}
\details{The function calculates Percent Intensity Difference (PID) 
or Percent Coefficient of Variation (CV) if there are two, or 
more than two replicates, respectively. PID is defined as percent ratio of the
absolute difference of replicate intensities and the mean of replicate intensities. 
CV is defined as percent ratio of standard deviation of 
replicate intensities and mean of replicate intensities.}
\value{Matrix with summary of feature consistency (min,first quartile (25th percentile), mean,median, 
third quartile (75th percentile), max,numgoodsamples, and QRscore).
QRscore=Percent good samples/median PID or CV}
\author{Karan Uppal <kuppal2@emory.edu>}
\keyword{ ~PID }
\keyword{ ~CV }
