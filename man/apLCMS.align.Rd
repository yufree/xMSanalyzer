\name{apLCMS.align}
\alias{apLCMS.align}
\title{apLCMS.align}
\description{Call apLCMS function, cdf.to.ftr, at different parameter settings}
\usage{
apLCMS.align(cdfloc, apLCMS.outloc, min.run.list = c(3,3), min.pres.list = c(0.3, 0.8), 
minexp = 2, mztol = 1e-5, alignmztol = NA, alignchrtol=NA, 
numnodes = 2, run.order.file = NA, subs = NA,filepattern=".cdf",
apLCMSmode = "untargeted", known_table, match_tol_ppm = 5)
}
\arguments{
  \item{cdfloc}{The folder where all NetCDF/mzML/mzXML/mzData files to be processed are located. For example "C:/experiment1/cdf/"}
  \item{apLCMS.outloc}{The folder where alignment output will be written. For example "C:/experiment1/apLCMSoutput/"}
  \item{min.run.list}{List of values for min.run parameter, eg: c(3,6) would run the cdf.to.ftr function at min.run=3 and min.run=6}
  \item{min.pres.list}{List of values min.pres, eg: c(0.3,0.8) would run the cdf.to.ftr function at min.run=3 and min.run=6}
  \item{minexp}{If a feature is to be included in the final feature table, it must be present in at least this number of spectra, eg:2}
   \item{mztol}{The user can provide the m/z tolerance level for peak identification to override the programs selection of the
   tolerance level. This value is expressed as the percentage of the m/z value. This value, multiplied by the m/z
   value, becomes the cutoff level. Please see the help for proc.cdf() for details.}
   \item{alignmztol}{The user can provide the m/z tolerance level for peak alignment to override the programs selection. This 
value is expressed as the percentage of the m/z value. This value, multiplied by the m/z value, becomes 
the cutoff level.Please see the help for feature.align() for details.}
  \item{alignchrtol}{The retention time tolerance level for peak alignment. 
  The default is NA, which allows the program to search for the tolerance level based on the data.
  Default: 10}
  \item{numnodes}{Number of nodes to use for processing.}
  \item{run.order.file}{Name of a tab-delimited file that includes sample names sorted by the order in which they were
  run (sample names must match the CDF file names)}
  \item{subs}{If not all the CDF files in the folder are to be processed, the user can define a subset using this parameter.
  For example, subs=15:30, or subs=c(2,4,6,8)}
   \item{filepattern}{File format of processed data files. eg: ".cdf", ".mzXML"}
    \item{apLCMSmode}{"untargeted" or "hybrid"; Default: "untargeted"}
  \item{known_table}{A data frame containing the known metabolite ions and previously found features. 
  Please see documentation of semi.sup() function in apLCMS for more details.}
  \item{match_tol_ppm}{The ppm tolerance to match identified features to known metabolites/features.
  Used by the semi.sup() function in apLCMS. Default: 5}
}
\details{This utility calls the cdf.to.ftr() function in the apLCMS package 
and performs serial sample processing at multiple combinations of two parameters:
min.run (minimum length of elution time for a series of signals grouped by m/z to be
considered a feature; default value: 3) and min.pres (minimum proportion of scans
in which the signal was present; default values: 0.3, 0.8). The function allows the
user to define parameters such as min.exp (minimum number of samples in which
a feature is present). This differs from the original apLCMS in that the original only
allows one set of parameters, whereas this function allows multiple sets. The
resulting tables containing m/z, retention time, and peak intensities in each
sample are stored at each parameter combination.}
\value{Feature table after weak signal recovery. This is the end product of the function cdf.to.ftr.}
\references{http://www.sph.emory.edu/apLCMS/ }
\author{Karan Uppal <kuppal2@emory.edu>}
\keyword{ ~apLCMS}
