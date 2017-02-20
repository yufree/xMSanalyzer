\name{find.Overlapping.mzs}
\alias{find.Overlapping.mzs}
\title{find.Overlapping.mzs}
\description{
This function matches features between two datasets using the
following user defined criteria:
1) Maximum m/z difference (+/-) ppm
2) Maximum retention time difference in seconds
}
\usage{find.Overlapping.mzs(dataA, dataB, mz.thresh = 10, time.thresh = NA,
alignment.tool=NA)}
\arguments{
  \item{dataA}{ apLCMS or XCMS feature table for dataset A}
  \item{dataB}{ apLCMS or XCMS feature table for dataset B}
  \item{mz.thresh}{Maximum m/z difference (+/-) ppm. eg: 10}
  \item{time.thresh}{Maximum retention time difference (+/-) secs. eg: 300}
  \item{alignment.tool}{Name of the feature alignment tool eg: "apLCMS" or "XCMS" or "NA"
  Use "NA" if the input matrix includes only m/z 
  or both m/z and retnetion time values.}
}
\details{
The find.Overlapping.mzs function operates on two sets of feature lists with m/z
and retention times for each feature, denoted by L1 and L2, and iterates over all
m/z values in L1 to find those that are within a user defined m/z (ppm) and
retention time (sec) threshold in L2. Optionally, the user can match features
based on only the m/z values by setting time.thresh=NA. The find.Unique.mzs
function uses a similar algorithm to find unique features that are
not within a user defined mass and retention time threshold level.
}
\value{
Matrix of overlapping features with columns:
index.data.A: index of overlapping m/z in dataset A
mz.data.A: m/z in dataset A
time.data.A: retention time in dataset A
index.data.B: index of overlapping m/z in dataset B
mz.data.B: m/z in dataset B
time.data.B: retention time in dataset B
}
\author{Karan Uppal <kuppal2@emory.edu>}
\seealso{apLCMS.align, XCMS.align, find.Unique.mzs, getVenn}
\keyword{ ~overlap }
\keyword{ ~common }