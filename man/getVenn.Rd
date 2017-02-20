\name{getVenn}
\alias{getVenn}
\title{getVenn}
\description{This utility calls the find.Overlapping.mzs function and generates a 
Venn diagram showing the extent of overlap between two datasets.}
\usage{getVenn(dataA, name_a, dataB, name_b, mz.thresh = 10, time.thresh=30,
alignment.tool, xMSanalyzer.outloc,use.unique.mz=FALSE,plotvenn=TRUE)}
\arguments{
  \item{dataA}{apLCMS or XCMS feature table for dataset A as a data frame.}
  \item{name_a}{Name of dataset A (eg: "SetA"). }
  \item{dataB}{apLCMS or XCMS feature table for dataset B as a data frame.}
  \item{name_b}{Name of dataset A (eg: "SetB")}
  \item{mz.thresh}{+/- ppm mass tolerance for m/z matching}
   \item{time.thresh}{Maximum retention time difference (+/-) secs. eg: 30}
     \item{alignment.tool}{"apLCMS" or "XCMS" or "NA". If NA is specified
     then the first two columns are treated as m/z and retention time,
     respectively.}
  \item{xMSanalyzer.outloc}{Output folder, eg: "C:/experiment1/xMSanalyzeroutput/"}
  \item{use.unique.mz}{If "TRUE", the function first finds unique features within each
  set}
   \item{plotvenn}{If "TRUE", the function plots the venn diagram.}  
}
\details{
This utility calls the find.Overlapping.mzs and find.Unique.mzs functions 
and generates a Venn diagram showing the extent of overlap between
datasets (up to three).
}
\value{A list is returned.
  \item{common }{Row numbers, m/zs, delta retention times of features 
  that are common between the two datasets.}
  \item{commonA }{Overlapping features in dataset A}
  \item{uniqueA }{Features that are unique in dataset A}
  \item{commonB }{Overlapping features in dataset B}
  \item{uniqueB }{Features that are unique in dataset B}
  \item{vennCounts}{Output of vennCounts function in limma package}
}
\references{http://rss.acs.unt.edu/Rdoc/library/limma/html/venn.html}
\author{Karan Uppal <kuppal2@emory.edu>}
\note{Only the unqiue m/zs within each dataset are used to generate Venn diagram.}
\keyword{ ~Venn }