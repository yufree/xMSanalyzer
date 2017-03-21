#' evaluate.Samples
#' 
#' Evaluate sample consistency based on Pearson or Spearman Correlation.
#' 
#' If at least two analytical replicates are present for each biological
#' sample, this function calculates the mean pairwise Pearson correlation
#' coefficient between sample replicates using the built-in cor() function in
#' R. Only the features with no missing values are used to evaluate
#' correlation. Analytical replicates refer to multiple injections from the
#' same biological sample; whereas, samples refer to different biological
#' samples.
#' 
#' @param curdata feature alignment output matrix from apLCMS or XCMS with
#' intensities
#' @param numreplicates number of technical replicates per sample
#' @param alignment.tool name of the feature alignment tool eg: "apLCMS" or
#' "XCMS"
#' @param cormethod Pearson or Spearman correlation.
#' @param missingvalue How are missing values represented? eg: 0 or NA
#' @param ignore.missing Should the missing values be ignored while computing
#' pearson correlation? eg: TRUE or FALSE
#' @param replace.bad.replicates Should the bad replicates be replaced by the
#' average of the good ones? For example, if the number of technical replicates
#' is more than two, and one of the replicates is poorly correlated with the
#' other two but the other replicates have correlation greater than the defined
#' threshold, then the bad replicate is replaced by the average of the good
#' ones.
#' @return returns a matrix of Pearson or Spearman Correlation Coefficients
#' within technical replicates per sample.
#' @author Karan Uppal <kuppal2@@emory.edu>
#' @keywords ~Sample quality
evaluate.Samples <- function(curdata, numreplicates, 
    alignment.tool, cormethod = "pearson", missingvalue = 0, 
    ignore.missing = TRUE, replace.bad.replicates = TRUE) {
    
    if (alignment.tool == "apLCMS") {
        col_end = 4
    } else {
        if (alignment.tool == "XCMS") {
            col_end = 8
        } else {
            stop(paste("Invalid value for alignment.tool. Please use either \"apLCMS\" or \"XCMS\"", 
                sep = ""))
        }
    }
    
    rnames <- colnames(curdata)
    rnames <- rnames[-c(1:col_end)]
    rnames = rnames[seq(1, length(rnames), numreplicates)]
    finalmat = {
    }
    
    curdata_int = curdata[, -c(1:col_end)]
    numsamp = dim(curdata_int)[2]
    
    curdata <- apply(curdata, 2, as.numeric)
    
    for (samp in seq(1, numsamp, numreplicates)) {
        samp_rep_last = samp + numreplicates - 1
        subdata = curdata_int[, samp:samp_rep_last]
        # if(FALSE)
        if (ignore.missing == TRUE) {
            if (is.na(missingvalue) == FALSE) {
                
                subdata <- replace(as.matrix(subdata), 
                  which(subdata == missingvalue), 
                  NA)
                
            }
            
        }
        rmat = cor(subdata, method = cormethod, use = "pairwise.complete.obs")
        
        rmat_upper = rmat[upper.tri(rmat)]
        
        good_reppairs <- which(rmat_upper > 0.7)
        check_bad_rep <- length(good_reppairs)
        
        rmat2 <- rmat
        diag(rmat2) <- 0
        bad_reppair <- apply(rmat2, 2, max)
        num_bad_reps <- length(which(bad_reppair < 
            0.7))
        
        if (numreplicates > 2) {
            if (num_bad_reps == 1) {
                
                bad_rep <- which(bad_reppair == min(bad_reppair))
                
                subdata[, c(bad_rep)] <- apply(subdata[, 
                  -c(bad_rep)], 1, function(x) {
                  mean(x, na.rm = TRUE)
                })
                
            }
            
        }
        
        if (replace.bad.replicates == TRUE) {
            curdata_int[, samp:samp_rep_last] <- subdata
            rmat = cor(subdata, method = cormethod, 
                use = "pairwise.complete.obs")
            
            rmat_upper = rmat[upper.tri(rmat)]
        }
        if (numreplicates == 2) {
            finalmat <- rbind(finalmat, mean(rmat_upper, 
                na.rm = TRUE))
        } else {
            if (numreplicates > 2) {
                rmat_vec = c(rmat_upper, mean(rmat_upper, 
                  na.rm = TRUE))
                finalmat <- rbind(finalmat, rmat_vec)
            }
            
        }
        
    }
    if (numreplicates == 2) {
        colnames(finalmat) <- c(paste(cormethod, "Correlation", 
            sep = ""))
        cnames <- "PearsonCorrelation"
    } else {
        if (numreplicates > 2) {
            
            cnames = {
            }
            for (repnum in seq(1, numreplicates - 
                1, 1)) {
                for (r1 in seq(repnum + 1, numreplicates, 
                  1)) {
                  cnames <- c(cnames, paste("rep", 
                    repnum, "vs", "rep", r1, sep = ""))
                }
            }
            cnames <- c(cnames, paste("mean", "Correlation", 
                sep = ""))
        }
    }
    if (replace.bad.replicates == TRUE) {
        curdata <- cbind(curdata[, c(1:col_end)], 
            curdata_int)
        curdata <- replace(as.matrix(curdata), which(is.na(curdata) == 
            TRUE), 0)
    }
    
    colnames(finalmat) <- cnames
    rownames(finalmat) <- rnames
    return(list(cor.matrix = finalmat, feature.table = curdata))
}
