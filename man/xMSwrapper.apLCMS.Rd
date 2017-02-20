\name{xMSwrapper.apLCMS}
\alias{xMSwrapper.apLCMS}
\title{xMSwrapper.apLCMS}
\description{Wrapper function based on apLCMS.align,evaluate.Features,
evaluate.Samples,merge.Results, and feat.batch.annotation.}
\usage{
xMSwrapper.apLCMS(cdfloc, apLCMS.outloc, xMSanalyzer.outloc,
                min.run.list = c(4,3), min.pres.list = c(0.5, 0.8),
                minexp.pct = 0.1, mztol = NA, alignmztol = NA,
                alignchrtol = NA, numnodes = NA, run.order.file = NA,
                apLCMSmode = "untargeted", known_table = NA,
                match_tol_ppm = 5, max.mz.diff = 15, max.rt.diff =
                300, merge.eval.pvalue = 0.2, mergecorthresh = 0.7,
                deltamzminmax.tol = 10, subs = NA, num_replicates = 3,
                mz.tolerance.dbmatch = 10, adduct.list = c("M+H"),
                samp.filt.thresh = 0.7, feat.filt.thresh = 50,
                cormethod = "pearson", mult.test.cor = TRUE,
                missingvalue = 0, ignore.missing = TRUE, filepattern =
                ".cdf", sample_info_file = NA, refMZ = NA,
                refMZ.mz.diff = 10, refMZ.time.diff = NA,
                void.vol.timethresh = 30, replacezeroswithNA = TRUE,
                scoreweight = 30,charge_type="pos", syssleep=0.5)
}
\arguments{
  \item{cdfloc}{The folder where all NetCDF/mzML/mzXML/mzData) files to be processed are located. For example "C:/CDF/"}
  \item{apLCMS.outloc}{The folder where apLCMS feature alignment output will be written. For example "C:/apLCMSoutput/"}
  \item{xMSanalyzer.outloc}{The folder where xMSanalyzer output will be written. For example "C:/xMSanalyzeroutput/"}
  \item{min.run.list}{List of values for min.run parameter, eg: c(3,6) would run the cdf.to.ftr function at 
  min.run=3 and min.run=6}
  \item{min.pres.list}{List of values min.pres, eg: c(0.3,0.8) would run the cdf.to.ftr function at min.run=3 and min.run=6}
  \item{minexp.pct}{If a feature is to be included in the final feature table, it must be present in at least this fraction of spectra. eg: 0.2}
  \item{mztol}{The user can provide the m/z tolerance level for peak identification to override the programs selection of the
tolerance level. This value is expressed as the percentage of the m/z value. This value, multiplied by the m/z
value, becomes the cutoff level. Please see the help for proc.cdf() for details.}
  \item{alignmztol}{The user can provide the m/z tolerance level for peak alignment to override the programs selection. This 
value is expressed as the percentage of the m/z value. This value, multiplied by the m/z value, becomes 
the cutoff level.Please see the help for feature.align() for details.}
  \item{numnodes}{Number of computational nodes to use for sample processing. eg: 2}
  \item{run.order.file}{Name of a tab-delimited file that includes sample names sorted by the order in which
  they were run (sample names must match the CDF file names)}
  \item{max.mz.diff}{+/- mz tolerance in ppm for feature matching}
  \item{max.rt.diff}{retention time tolerance for feature matching}
  \item{merge.eval.pvalue}{Threshold for defining signifcance level of the paired t-test or
  the Pearson correlation during the merge stage in xMSanalyzer. The p-value is used
  to determine whether two features with same m/z and retention time have identical
  intensity profiles.}
  \item{mergecorthresh}{Correlation threshold to be used during the merge stage in xMSanalyzer
  to determine whether two features with same m/z and retention time have identical
  intensity profiles.}
  \item{subs}{If not all the CDF files in the folder are to be processed, the user can define a subset using this parameter. 
 For example, subs=15:30, or subs=c(2,4,6,8)}
  \item{num_replicates}{Number of replicates per sample}
  \item{mz.tolerance.dbmatch}{m/z threshold for database matching}
  \item{adduct.list}{List of adducts for matching m/zs in KEGG. eg: c("M+H", "M+Na")}
  \item{samp.filt.thresh}{Threshold for filtering samples based on Pearson correlation between technical replicates. eg: 0.7}
  \item{feat.filt.thresh}{Threshold for filtering samples based on percent intensity difference or coefficient of variation. eg: 50}
  \item{cormethod}{Method for determing correlation between technical
  replicates. eg: "pearson" or "spearman}
   \item{mult.test.cor}{Should Bonferroni multiple testing correction method 
  be applied for comparing intensities of overlapping m/z? Default: FALSE}
  \item{missingvalue}{How are missing values represented? eg: 0 or NA}
  \item{ignore.missing}{Should the missing values be ignored while computing
  pearson correlation. eg: TRUE or FALSE}
  \item{filepattern}{File format of spectral data files. eg: ".cdf", ".mzXML"}
  \item{alignchrtol}{The retention time tolerance level for peak alignment. 
  The default is NA, which allows the program to search for the tolerance level based on the data.
  Default: 10}
  \item{apLCMSmode}{"untargeted" or "hybrid"; Default: "untargeted"}
  \item{known_table}{A data frame containing the known metabolite ions and previously found features. 
  Please see documentation of semi.sup() function in apLCMS for more details.}
  \item{match_tol_ppm}{The ppm tolerance to match identified features to known metabolites/features.
  Used by the semi.sup() function in apLCMS. Default: 5}
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
  \item{void.vol.timethresh}{Threshold for void volume. The program searches for the void volume
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
  \item{charge_type}{Ionization mode. "pos" or "neg" Default: "pos"}
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
  \item{apLCMS.merged.res }{Merged feature table, P1 U P2
 where P1 and P2 are two sets of parameter settings. Please note that the four columns
 include the characteristics of the feature from apLCMS. The "npeaks" column includes
 the number of samples with a non-missing intensity. The "QRscore" is defined as
 the ratio of percentage of biological samples for which more than 50 percent of technical replicates have
 a signal and median PID or CV. QRscore=Percent good samples/median PID or CV}
  \item{apLCMS.ind.res}{List with results from individual parameter settings.}
  \item{apLCMS.ind.res.filtered}{List with results from individual parameter settings after filtering.}
  \item{annot.res}{List with annotation results after merging results.}
  \item{feat.eval.ind}{List with feature evaluation results based on CV or PID
  at each parameter setting.}
  \item{sample.eval.ind}{List with sample evaluation results at each parameter setting based on
  correlation between technical replicates.}
Outputs the following to apLCMS.outloc:
-apLCMS results at each parameter settings
Outputs the following to xMSanalyzer.outloc:
-feature consistency results
-sample evaluation results
-feature list after merging results from parameter settings A and B
-merge summary
}
\author{Karan Uppal <kuppal2@emory.edu>}
\keyword{ ~apLCMS }
\keyword{ ~xMSwrapper }
