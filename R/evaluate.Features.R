#' evaluate.Features
#' 
#' Evaluates feature consistency within analytical/technical replicates of
#' samples based on PID or CV. Reports quantile distribution of PID or CV
#' within technical replicates across all samples and computes a quantitative
#' reproducibility score, QRscore, which is defined as the ratio of percentage
#' of biological samples for which more than 50 percent of technical replicates
#' have a signal and median PID or CV
#' 
#' The function calculates Percent Intensity Difference (PID) or Percent
#' Coefficient of Variation (CV) if there are two, or more than two replicates,
#' respectively. PID is defined as percent ratio of the absolute difference of
#' replicate intensities and the mean of replicate intensities.  CV is defined
#' as percent ratio of standard deviation of replicate intensities and mean of
#' replicate intensities.
#' 
#' @param curdata Feature table from apLCMS or XCMS with sample intensities.
#' @param numreplicates Number of replicates per sample. eg: 2
#' @param min.samp.percent If signal is detected in x proportion of technical
#' replicates, then the missing values are replaced by mean intensity which is
#' calculated using non-missing replicates.  Defaul: 0.6
#' @param alignment.tool Name of the feature alignment tool eg: "apLCMS" or
#' "XCMS".
#' @param impute.bool Logical (TRUE or FALSE). Used for calculating coefficient
#' of variation (CV). Missing values are replaced by the mean of the
#' non-missing intensities if at least 60% of analytical replicates have
#' values. For example, if two out of three replicates have non-missing values
#' (intensity >0) and the third replicate has a missing value (intensity=0),
#' then the intensity of the third replicate is replaced by the mean of the
#' intensities of the other two to calculate CV.
#' @param missingvalue How are the missing values represented? eg: 0 or NA
#' @return Matrix with summary of feature consistency (min,first quartile (25th
#' percentile), mean,median, third quartile (75th percentile),
#' max,numgoodsamples, and QRscore). QRscore=Percent good samples/median PID or
#' CV
#' @author Karan Uppal <kuppal2@@emory.edu>
#' @keywords ~PID ~CV
evaluate.Features <- function(curdata, numreplicates, 
    min.samp.percent = 0.6, alignment.tool, impute.bool, 
    missingvalue = 0) {
    
    if (numreplicates == 2) {
        print("**calculating percent intensity difference**")
        eval.feat.results <- getPID(curdata, alignment.tool, 
            missingvalue)
    } else {
        if (numreplicates > 2) {
            print("**calculating CV**")
            eval.feat.results <- getCVreplicates(curdata, 
                alignment.tool, min.samp.percent = 0.6, 
                numreplicates, impute.bool, missingvalue)
        }
        
    }
    return(eval.feat.results)
}
