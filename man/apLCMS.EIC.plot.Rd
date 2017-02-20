\name{apLCMS.EIC.plot}
\alias{apLCMS.EIC.plot}
\title{apLCMS.EIC.plot}
\description{Modified version of the EIC plot function in apLCMS}
\usage{
apLCMS.EIC.plot(aligned, rows = NA, colors = NA, transform = "none", subset = NA, 
minrt = NA, maxrt = NA, min.run = 12, min.pres = 0.5, max.spline.time.points = 1000)
}
\arguments{
  \item{aligned}{Rda object from apLCMS}
  \item{rows}{feature row indices eg: c(10,35)}
  \item{colors}{color vector eg: c("red","green")}
  \item{transform}{Tranformation method eg: "none", "log", "sqrt", or "cuberoot"}
  \item{subset}{subset of profiles for which to plot the EICs eg: c(1:6)}
  \item{minrt}{minimum retention limit in EIC plot eg:30}
  \item{maxrt}{maximum retention time limit in EIC plot eg: 300}
  \item{min.run}{same as apLCMS.align}
  \item{min.pres}{same as apLCMS.align}
  \item{max.spline.time.points}{Maximum number of time points to use for
	interpolation spline eg: 1000}
}
\details{
This function is modified version of the EIC.plot function in apLCMS.
Users can define the time range for plotting EICs.
}
\author{Karan Uppal <kuppal2@emory.edu>; Tianwei Yu}
\keyword{ ~apLCMS }
\keyword{ ~EIC }