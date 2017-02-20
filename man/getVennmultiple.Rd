\name{getVennmultiple}
\alias{getVennmultiple}
\title{getVennmultiple}
\description{This utility calls the find.Overlapping.mzs function and generates a
Venn diagram showing the extent of overlap between three datasets.
}
\usage{
getVennmultiple(dataA, name_a, dataB, name_b, dataC, name_c, mz.thresh = 10, 
time.thresh = 30, alignment.tool, xMSanalyzer.outloc,use.unique.mz=FALSE,
plotvenn=TRUE)
}
\arguments{
  \item{dataA}{apLCMS or XCMS feature table for dataset A as a data frame.}
  \item{name_a}{Name of dataset A (eg: "SetA").}
  \item{dataB}{apLCMS or XCMS feature table for dataset B as a data frame.}
  \item{name_b}{Name of dataset B (eg: "SetB")}
  \item{dataC}{apLCMS or XCMS feature table for dataset B as a data frame.}
  \item{name_c}{Name of dataset C (eg: "SetC")}
  \item{mz.thresh}{+/- ppm threshold for m/z matching. eg: 10}
  \item{time.thresh}{Maximum retention time difference (+/-) secs. eg: 30}
 \item{alignment.tool}{Name of the feature alignment tool eg: "apLCMS" or "XCMS" or "NA"
  Use "NA" if the input matrix includes only m/z 
  or both m/z and retnetion time values.}
  \item{xMSanalyzer.outloc}{xMSanalyzer output location, 
  eg: "C:/experiment1/xMSanalyzeroutput/}
   \item{use.unique.mz}{If "TRUE", the function first finds unique features within each
  set}
   \item{plotvenn}{If "TRUE", the function plots the venn diagram.}  
}
\value{A list is returned.
 \item{uniqueA }{Features that are unique in dataset A}
  \item{uniqueB }{Features that are unique in dataset B}
  \item{uniqueC }{Features that are unique in dataset C}
  \item{common }{Features that are common}
  \item{vennCounts}{Output of vennCounts function in limma package}
}
\references{http://rss.acs.unt.edu/Rdoc/library/limma/html/venn.html }
\author{Karan Uppal <kuppal2@emory.edu>}
\note{Only the unqiue m/zs within each dataset are used to generate Venn diagram.}
\keyword{ ~Venn }
