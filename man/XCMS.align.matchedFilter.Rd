\name{XCMS.align.matchedFilter}
\alias{XCMS.align.matchedFilter}
\title{XCMS.align.matchedFilter}
\description{Runs XCMS using the matchedFilter alignment algorithm at different parameter settings.}
\usage{
XCMS.align.matchedFilter(cdfloc, XCMS.outloc, step.list = c(0.001), mz.diff.list = c(0.1), 
sn.thresh.list = c(3), max =50, bw.val = c(10), minfrac.val = 0.5, minsamp.val = 2, 
mzwid.val = 0.25, sleep.val = 0, run.order.file = NA, subs = NA, retcor.family = "symmetric",
retcor.plottype = "mdevden", groupval.method = "medret")
}
\arguments{
\item{cdfloc}{The folder where all CDF/mzXML files to be processed are located. For example "C:/CDF/"}
  \item{XCMS.outloc}{The folder where alignment output will be written. For example "C:/CDFoutput/"}
  \item{step.list}{list containing values for the step size}
  \item{mz.diff.list}{list containing values for the minimum difference for features with retention time overlap}
  \item{sn.thresh.list}{list containing values for signal to noise ratio cutoff variable}
  \item{max}{Value for maxnimum number of peaks per EIC variable: eg: 5}
  \item{bw.val}{bandwidth value}
  \item{minfrac.val}{minimum fraction of samples necessary in at least one of the sample groups for it to be
  a valid group}
  \item{minsamp.val}{minimum number of samples necessary in at least one of the sample groups for it to be
  a valid group}
  \item{mzwid.val}{width of overlapping m/z slices to use for creating peak density chromatograms and
  grouping peaks across samples}
  \item{sleep.val}{seconds to pause between plotting successive steps of the
peak grouping algorithm. peaks are plotted as points showing
relative intensity. identified groups are flanked by dotted
vertical lines.}
  \item{run.order.file}{Name of a tab-delimited file that includes sample names sorted by the order in
  which they were run(sample names must match the CDF file names)}
  \item{subs}{If not all the CDF files in the folder are to be processed, the user can define a subset
  using this parameter. For example, subs=15:30, or subs=c(2,4,6,8)}
 \item{retcor.family}{Used by matchedFilter alignment method. Use "gaussian" to perform fitting by
				least squares without outlier removal. Or "symmetric" to use a redescending
				M estimator with Tukey's biweight function that allows outlier removal.}
  \item{retcor.plottype}{Used by both matchedFilter and centWave alignment methods. eg: "deviation"
				or "mdevden"}
  \item{groupval.method}{Conflict resolution method while calculating peak values for 
  each group. eg: "medret" or "maxint"}
}
\details{This is a wrapper function based on the xcms Bioconductor package for preprocessing/analysis
of mass spectral data. The XCMS.align utility performs serial sample processing at multiple combinations 
of four parameters: step (the step size; default values: 0.001, 0.01, 0.1),
mzdiff (minimum difference for features with retention time overlap; default
values: 0.001, 0.01, 0.1), snthresh (signal-to-noise ratio cutoff; default values:
3, 6, 10), and max (maximum number of peaks per EIC; default values: 5, 10).
The resulting tables containing m/z, retention time, and mean peak intensities
in each sample are stored at each parameter combination.}
\value{A matrix, with columns of m/z values, elution times, mean signal strengths in each spectrum}
\author{Karan Uppal <kuppal2@emory.edu>}
\note{Please refer to the xcms manual in Bioconductor for more details.}
\keyword{ ~alignment }
\keyword{ ~xcms }
\keyword{ ~matchedFilter}
