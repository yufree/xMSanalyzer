\name{XCMS.align.centWave}
\alias{XCMS.align.centWave}
\title{XCMS.align.centWave}
\description{Wrapper function for XCMS using the centwave alignment algorithm.}
\usage{
XCMS.align.centWave(cdfloc, XCMS.outloc, ppm.list = c(10), mz.diff.list = c(-0.001), 
sn.thresh.list = c(10), prefilter.list = c(3, 100), bw.val = c(10), 
groupval.method = "medret", step.list = c(0.1), max = 50, minfrac.val = 0.5, 
minsamp.val = 1, mzwid.val = 0.25, sleep.val = 0, run.order.file = NA, subs = NA, 
retcor.method = "obiwarp", retcor.family = "symmetric", retcor.plottype = "deviation",
peakwidth = c(20, 50), nSlaves=2)
}
\arguments{
  \item{cdfloc}{The folder where all CDF/mzXML files to be processed are located. For example "C:/CDF/"}
  \item{XCMS.outloc}{The folder where alignment output will be written. For example "C:/CDFoutput/"}
  \item{ppm.list}{list containing values for maximal tolerated m/z deviation in consecutive scans, in ppm}
  \item{mz.diff.list}{list containing values for the minimum difference for features with
  retention time overlap. eg: c(0.001,0.1)}
  \item{sn.thresh.list}{list containing values for signal to noise ratio cutoff variable. eg: c(3,10)}
  \item{prefilter.list}{prefiltering values c(k,l) where mass traces that do not contain at least k peaks
  with intensity>=l are filtered}
   \item{bw.val}{bandwidth value}
  \item{groupval.method}{Conflict resolution method while calculating peak values for 
  each group. eg: "medret" or "maxint"}
 \item{step.list}{list containing values for the step size. eg: c(0.1,1)}
   \item{max}{Value for maxnimum number of peaks per EIC variable. eg: 50}
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
  \item{retcor.method}{Method for aligning retention times across samples. eg: "loess" or "obiwarp"}
  \item{retcor.family}{Used by matchedFilter alignment method. Use "gaussian" to perform fitting by
				least squares without outlier removal. Or "symmetric" to use a redescending
				M estimator with Tukey's biweight function that allows outlier removal.}
  \item{retcor.plottype}{Used by both matchedFilter and centWave alignment methods. eg: "deviation"
				or "mdevden"}
  \item{peakwidth}{Chromtagrophic peak width in seconds. eg: c(20,50)}
  \item{nSlaves}{Number of computing cores to be used. eg: 2}
}
\details{This is a wrapper function based on the xcms Bioconductor package for preprocessing/analysis
of mass spectral data. The resulting tables containing m/z, retention time, and mean peak intensities
in each sample are stored at each parameter combination.}
\value{A matrix, with columns of m/z values, elution times, mean signal strengths in each spectrum}
\references{Tautenhahn R, Bottcher C, Neumann S. Highly sensitive feature detection for high resolution LC/MS. 
BMC Bioinformatics. 2008 Nov 28.}
\author{Karan Uppal}
\note{Please refer to the xcms manual in Bioconductor for more details.}
\keyword{ ~alignment }
\keyword{ ~xcms }
