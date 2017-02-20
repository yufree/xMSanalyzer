\name{xMSwrapper.XCMS.centWave}
\alias{xMSwrapper.XCMS.centWave}
\title{xMSwrapper.XCMS.centWave}
\description{Wrapper function for XCMS using the centwave alignment algorithm,
evaluate.Features,evaluate.Samples,merge.Results, and feat.batch.annotation}
\usage{
xMSwrapper.XCMS.centWave(cdfloc, XCMS.outloc, xMSanalyzer.outloc, 
ppm.list =c(10, 25, 30), mz.diff.list = c(-0.001, 0.1), sn.thresh.list = c(3, 5, 10),
prefilter.list = c(3,100), bw.val = c(10, 30), groupval.method = "medret",
step.list = c(0.1, 1), max = 50, minfrac.val = 0.5,minsamp.val = 2, mzwid.val = 0.25, sleep.val = 0,
retcor.method = "obiwarp", retcor.family ="symmetric", retcor.plottype = "deviation",
peakwidth = c(20, 50), numnodes = 2, run.order.file = NA,max.mz.diff = 15, 
max.rt.diff = 300, merge.eval.pvalue = 0.2, mergecorthresh = 0.7, 
deltamzminmax.tol = 10,num_replicates = 2, subs = NA, mz.tolerance.dbmatch =10,
adduct.list = c("M+H"), samp.filt.thresh = 0.7,feat.filt.thresh = 50,
cormethod = "pearson",mult.test.cor = TRUE, missingvalue = 0, ignore.missing= TRUE,
sample_info_file = NA, refMZ = NA,refMZ.mz.diff = 10, refMZ.time.diff = NA,
void.vol.timethresh = 30, replacezeroswithNA = TRUE, scoreweight = 30,
filepattern = ".cdf", charge_type="pos",minexp.pct=0.1,syssleep=0.5)
}
\arguments{
  \item{cdfloc}{The folder where all CDF/mzXML files to be processed are located. For example "C:/CDF/"}
  \item{XCMS.outloc}{The folder where alignment output will be written. For example "C:/CDFoutput/"}
  \item{xMSanalyzer.outloc}{The folder where xMSanalyzer output will be written. For example "C:/xMSanalyzeroutput/"}
   \item{ppm.list}{list containing values for maximal tolerated m/z deviation in consecutive scans, in ppm}
  \item{mz.diff.list}{list containing values for the minimum difference for features with retention time overlap}
  \item{sn.thresh.list}{list containing values for signal to noise ratio cutoff variable}
  \item{prefilter.list}{prefiltering values c(k,l) where mass traces that do not contain at least k peaks
  with intensity>=l are filtered}
   \item{bw.val}{bandwidth value}
  \item{groupval.method}{Conflict resolution method while calculating peak values for 
  each group. eg: "medret" or "maxint"}
 \item{step.list}{list containing values for the step size}
   \item{max}{list containing values for maxnimum number of peaks per EIC variable}
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
  \item{retcor.method}{Method for aligning retention times across samples. eg: "loess" or "obiwarp"}
  \item{retcor.family}{Used by matchedFilter alignment method. Use "gaussian" to perform fitting by
				least squares without outlier removal. Or "symmetric" to use a redescending
				M estimator with Tukey's biweight function that allows outlier removal.}
  \item{retcor.plottype}{Used by both matchedFilter and centWave alignment methods. eg: "deviation"
				or "mdevden"}
  \item{peakwidth}{Chromtagrophic peak width in seconds. eg: c(20,50)}
  \item{numnodes}{Number of computing nodes to use.  eg: 1}
  \item{run.order.file}{Name of a tab-delimited file that includes sample names sorted by the order in which they were run (sample names must match the CDF file names)}
  \item{max.mz.diff}{+/- mz tolerance in ppm for feature matching}
  \item{max.rt.diff}{retention time tolerance for feature matching}
  \item{merge.eval.pvalue}{Threshold for defining signifcance level of the paired t-test or
  the Pearson correlation during the merge stage in xMSanalyzer. The p-value is used
  to determine whether two features with same m/z and retention time have identical
  intensity profiles.}
  \item{mergecorthresh}{Correlation threshold to be used during the merge stage in xMSanalyzer
  to determine whether two features with same m/z and retention time have identical
  intensity profiles.}
  \item{num_replicates}{Number of replicates per sample}
  \item{subs}{If not all the CDF files in the folder are to be processed, the user can define a subset using this parameter.}
  \item{mz.tolerance.dbmatch}{m/z threshold for database matching}
  \item{adduct.list}{List of adducts for matching m/zs in KEGG. eg: c("M+H","M+H-H2O")}
  \item{samp.filt.thresh}{Threshold for filtering samples based on Pearson correlation between technical replicates. eg: 0.7}
  \item{feat.filt.thresh}{Threshold for filtering samples based on percent intensity difference or coefficient of variation. eg: 50}
  \item{cormethod}{Method for determing correlation between technical
  replicates. eg: "pearson" or "spearman}
  \item{mult.test.cor}{Should Bonferroni multiple testing correction method 
  be applied for comparing intensities of overlapping m/z? Default: FALSE}
  \item{missingvalue}{How are missing values represented? eg: 0 or NA}
  \item{ignore.missing}{Should the missing values be ignored while computing
  pearson correlation. eg: TRUE or FALSE}
   \item{deltamzminmax.tol}{Maximum allowed delta ppm between mz.min and mz.max. Eg: 10}
  \item{refMZ}{Full path of the file with m/z of the targeted chemicals to search for.
    If the value is "NA", then the list of metabolites in the data(example_target_list) is used.
    The input file should be formatted as data(example_target_list):
  Column A: m/z of the targeted metabolite (required)
  Column B: retention time (Optional)
  Column C: Name of the metabolite (required)
  }
  \item{refMZ.mz.diff}{The ppm tolerance to search for targeted metabolites/features in xMSanalyzer.
  Default: 10}
  \item{refMZ.time.diff}{The time tolerance to search for targeted metabolites/features in xMSanalyzer.
  Default: NA}
  \item{void.vol.timethresh}{Time threshold for void volume. The program searches for the void volume
  cutoff within the defined time limit. Default: 30}
  \item{replacezeroswithNA}{Should 0s be treated as missing values during ComBat (TRUE or FALSE).
  Default: TRUE} 
  \item{scoreweight}{The w parameter in the scoring function defined in the xMSanalyzer manuscript.
  Uppal 2013, BMC Bioinformatiocs. Default: 30}
  \item{sample_info_file}{File listing the order in which the samples were run. The format of the file
  should be as follows:
  Column A:File names matching the CDF/mzXML files
  Column B: Sample ID
  Column C: Batch number (This column should be labeled "Batch")
  Column D: Additional covariates to adjust for (Optional)
  }
  \item{filepattern}{File format of spectral data files.
  eg: ".cdf", ".mzXML"}
  \item{charge_type}{Ionization mode. "pos" or "neg" Default: "pos"}
  \item{minexp.pct}{Minimum fraction of samples in which the signal should be
  present. Default: "0.1"}
  \item{syssleep}{Sleep time in between KEGG REST queries. Default: "0.5"}
  
}
\details{The wrapper function includes five stages to utilize information from the technical replicates
to optimize the data extraction process, enhance data quality, search for targeted features/metabollites,
perform QA & QC evaluations including batch-effect evaluation & correction. :
1) features are extracted using different parameters 
2) results from each parameter setting are evaluated for sample quality & feature consistency, 
3) filtered results from individual settings are merged to obtain a combined feature table; the 
optimization score that takes into account the number of features and average CV (see Uppal 2013)
is used to determine the most optimal set of parameters.
4) A targeted feature table using the list of m/z in the "refMZ" file is generated
5) Quality measures of each featue includes: number of samples including replicates 
with non-missing values (NumPres.All.Samples), number of  biological samples for which more than
60% of technical replicates have a non-missing/non-zero value (NumPres.Biol.Samples),
median coefficient of variation (CV) within technical replicates summarized across all samples; 
Qscore (Quality score which is calculated using NumPres.All.Samples, NumPres.Biol.Samples,
median CV, and delta ppm between mz.min and mz.max; Higher is better) 
Users have the option to filter poor quality samples and features based on correlation between technical
replicates and feature reproducibility measures such as PID or CV, respectively.}
\value{A list is returned.
  \item{XCMS.merged.res }{Merged feature table, P1 U P2
 where P1 and P2 are two sets of parameter settings. The "QRscore" is defined as
 the ratio of percentage of biological samples for which more than 50 percent of technical replicates have
 a signal and median PID or CV. QRscore=Percent good samples/median PID or CV}
  \item{XCMS.ind.res}{List with results from individual parameter settings.}
  \item{XCMS.ind.res.filtered}{List with results from individual parameter settings after filtering.}
  \item{annot.res}{List with annotation results after merging results.}
  \item{feat.eval.ind}{List with feature evaluation results based on CV or PID
  at each parameter setting.}
  \item{sample.eval.ind}{List with sample evaluation results at each parameter setting based on
  correlation between technical replicates.}
Outputs XCMS results at each parameter settings to XCMS.outloc and the following to xMSanalyzer.outloc:
 -feature consistency results
-sample evaluation results
-feature list after merging results from parameter settings A and B
-merge summary
}
\author{Karan Uppal <kuppal2@emory.edu>}
\keyword{ ~XCMS }
\keyword{ ~centWave }
