\name{find.Unique.mzs}
\alias{find.Unique.mzs}
\title{find.Unique.mzs}
\description{This function finds unique m/zs between two datasets.}
\usage{find.Unique.mzs(dataA, dataB, mz.thresh = 10, time.thresh=NA, alignment.tool)}
\arguments{
  \item{dataA}{apLCMS or XCMS feature table for dataset A}
  \item{dataB}{apLCMS or XCMS feature table for dataset B}
  \item{mz.thresh}{Maximum m/z difference (+/-) ppm}
  \item{time.thresh}{Maximum retention time difference (+/-) seconds}
  \item{alignment.tool}{Name of the feature alignment tool eg: "apLCMS" or "XCMS" or "NA"
  Use "NA" if the input matrix includes only m/z 
  or both m/z and retnetion time values.}
}
\details{
The find.Unique.mzs function operates on two sets of feature lists with m/z for each feature,
denoted by L1 and L2, and iterates over all m/z values in L1 to find those that are within a
user defined m/z (ppm) threshold in L2. 
}
\value{A list is returned:
  \item{uniqueA }{Unique m/zs in dataA}
  \item{uniqueB }{Unique m/zs in dataB}
}

\author{Karan Uppal <kuppal2@emory.edu>}
\seealso{apLCMS.align, XCMS.align, find.Overlapping.mzs}
\keyword{ ~unique }