#' XCMS.align.matchedFilter
#' 
#' Runs XCMS using the matchedFilter alignment algorithm at different parameter
#' settings.
#' 
#' This is a wrapper function based on the xcms Bioconductor package for
#' preprocessing/analysis of mass spectral data. The XCMS.align utility
#' performs serial sample processing at multiple combinations of four
#' parameters: step (the step size; default values: 0.001, 0.01, 0.1), mzdiff
#' (minimum difference for features with retention time overlap; default
#' values: 0.001, 0.01, 0.1), snthresh (signal-to-noise ratio cutoff; default
#' values: 3, 6, 10), and max (maximum number of peaks per EIC; default values:
#' 5, 10). The resulting tables containing m/z, retention time, and mean peak
#' intensities in each sample are stored at each parameter combination.
#' 
#' @param cdfloc The folder where all CDF/mzXML files to be processed are
#' located. For example "C:/CDF/"
#' @param XCMS.outloc The folder where alignment output will be written. For
#' example "C:/CDFoutput/"
#' @param step.list list containing values for the step size
#' @param mz.diff.list list containing values for the minimum difference for
#' features with retention time overlap
#' @param sn.thresh.list list containing values for signal to noise ratio
#' cutoff variable
#' @param max Value for maxnimum number of peaks per EIC variable: eg: 5
#' @param bw.val bandwidth value
#' @param minfrac.val minimum fraction of samples necessary in at least one of
#' the sample groups for it to be a valid group
#' @param minsamp.val minimum number of samples necessary in at least one of
#' the sample groups for it to be a valid group
#' @param mzwid.val width of overlapping m/z slices to use for creating peak
#' density chromatograms and grouping peaks across samples
#' @param sleep.val seconds to pause between plotting successive steps of the
#' peak grouping algorithm. peaks are plotted as points showing relative
#' intensity. identified groups are flanked by dotted vertical lines.
#' @param run.order.file Name of a tab-delimited file that includes sample
#' names sorted by the order in which they were run(sample names must match the
#' CDF file names)
#' @param subs If not all the CDF files in the folder are to be processed, the
#' user can define a subset using this parameter. For example, subs=15:30, or
#' subs=c(2,4,6,8)
#' @param retcor.family Used by matchedFilter alignment method. Use "gaussian"
#' to perform fitting by least squares without outlier removal. Or "symmetric"
#' to use a redescending M estimator with Tukey's biweight function that allows
#' outlier removal.
#' @param retcor.plottype Used by both matchedFilter and centWave alignment
#' methods. eg: "deviation" or "mdevden"
#' @param groupval.method Conflict resolution method while calculating peak
#' values for each group. eg: "medret" or "maxint"
#' @return A matrix, with columns of m/z values, elution times, mean signal
#' strengths in each spectrum
#' @note Please refer to the xcms manual in Bioconductor for more details.
#' @author Karan Uppal <kuppal2@@emory.edu>
#' @keywords ~alignment ~xcms ~matchedFilter
XCMS.align.matchedFilter <- function(cdfloc, XCMS.outloc, 
    step.list = c(0.001), mz.diff.list = c(0.1), sn.thresh.list = c(3), 
    max = 50, bw.val = c(10), minfrac.val = 0.5, minsamp.val = 2, 
    mzwid.val = 0.25, sleep.val = 0, run.order.file = NA, 
    subs = NA, retcor.family = "symmetric", retcor.plottype = "mdevden", 
    groupval.method = "medret", target.mz.list = NA) {
    
    setwd(cdfloc)
    dir.create(XCMS.outloc, showWarnings = FALSE)
    cdf_files = list.files(cdfloc, ".cdf|.mzxml|mXML", 
        ignore.case = TRUE)
    aligned_data_list = new("list")
    pcount = 1
    
    if (is.na(subs[1]) == FALSE) {
        cdf_files = cdf_files[subs]
        numsamp = length(subs)
        
        
    } else {
        
        numsamp = length(cdf_files)
    }
    
    for (t in sn.thresh.list) {
        for (s in step.list) {
            for (m in mz.diff.list) {
                # for(maxl in max.list)
                {
                  
                  xset = xcmsSet(cdf_files, step = s, 
                    snthresh = t, mzdiff = m, max = max)
                  
                  
                  xset <- group(xset)
                  
                  xset2 <- retcor(xset, family = retcor.family, 
                    plottype = retcor.plottype)
                  
                  for (b in bw.val) {
                    ### Group peaks together across samples, set
                    ### bandwitdh, change important m/z parameters here
                    ### Syntax: group(object, bw = 30, minfrac = 0.5,
                    ### minsamp= 1, mzwid = 0.25, max = 5, sleep = 0)
                    xset2 <- group.density(xset2, 
                      bw = b, minfrac = minfrac.val, 
                      minsamp = minsamp.val, mzwid = mzwid.val, 
                      max = max, sleep = sleep.val)
                    # xset3=fillPeaks(xset2)
                    xset3 = suppressWarnings(fillPeaks(xset2))
                    
                    print("Getting sample intensities")
                    finalfeatmat = {
                    }
                    
                    if (dim(xset3@groups)[1] > 0) {
                      groupmat <- groups(xset3)
                      
                      
                      # group_intmat<-groupval(xset3, method =
                      # groupval.method, intensity = 'into')
                      
                      group_intmat <- groupval(xset3, 
                        groupval.method, "into")
                      finalfeatmat <- cbind(groupmat, 
                        group_intmat)
                      finalfeatmat <- as.data.frame(finalfeatmat, 
                        row.names = NULL)
                      
                    } else {
                      if (length(xset3@sampnames) == 
                        1) {
                        finalfeatmat <- xset3@peaks
                      } else {
                        stop("First argument must be a xcmsSet with group information or contain only one sample.")
                      }
                    }
                    
                    
                    fname = paste(XCMS.outloc, "/XCMS_matchedFilter", 
                      "_thresh", t, "_step", s, "_mzdiff", 
                      m, "_max", max, "_bw", b, ".txt", 
                      sep = "")
                    cnames = colnames(finalfeatmat)
                    cnames[1] = "mz"
                    cnames[4] = "time"
                    colnames(finalfeatmat) = c(cnames[1:8], 
                      cdf_files)
                    cnames <- tolower(cnames)
                    cnames <- gsub(".cdf", "", cnames)
                    
                    if (is.na(run.order.file) == FALSE) {
                      fileorder = read.table(run.order.file, 
                        header = FALSE)
                      fileorder = apply(fileorder, 
                        1, tolower)
                      ordlist = sapply(1:length(fileorder), 
                        function(i) {
                          which(cnames == fileorder[i])
                        })
                      ordlist = unlist(ordlist)
                      ordlist = ordlist + 8
                      
                      finalfeatmat = finalfeatmat[, 
                        c(1:8, ordlist)]
                    }
                    
                    write.table(finalfeatmat, file = fname, 
                      sep = "\t", row.names = FALSE)
                    aligned_data_list[[pcount]] <- finalfeatmat
                    pcount = pcount + 1
                    
                    if (is.na(target.mz.list) == FALSE) {
                      
                      fname = paste(XCMS.outloc, "/EICmatchedFilter", 
                        "_thresh", t, "_step", s, 
                        "_mzdiff", m, "_max", max, 
                        "_bw", b, ".pdf", sep = "")
                      pdf(fname)
                      if (is.na(target.mz.list[1, 
                        1]) == FALSE) {
                        stddata <- target.mz.list
                        
                        # print(head(stddata))
                      } else {
                        Name <- paste("mz", seq(1, 
                          dim(finalfeatmat)[1]), sep = "")
                        stddata <- cbind(finalfeatmat[, 
                          c(1)], Name)
                        
                      }
                      overlapres5ppm <- getVenn(dataA = finalfeatmat, 
                        name_a = paste("Expdata", 
                          sep = ""), dataB = stddata, 
                        name_b = "Target", mz.thresh = 10, 
                        time.thresh = NA, alignment.tool = "XCMS", 
                        xMSanalyzer.outloc = XCMS.outloc, 
                        plotvenn = FALSE)
                      if (length(unique(overlapres5ppm$common$index.A)) > 
                        0) {
                        setwd(cdfloc)
                        num_samples <- dim(finalfeatmat)[2] - 
                          8
                        min_samp <- 6
                        if (num_samples < min_samp) {
                          min_samp <- num_samples
                        }
                        rand_sample_set <- sample(size = num_samples, 
                          x = num_samples, replace = FALSE)
                        rand_set <- c(1:min_samp, 
                          (num_samples - 2):num_samples)
                        com_ind <- which(rand_sample_set %in% 
                          rand_set)
                        if (length(com_ind) > 0) {
                          rand_sample_set <- rand_sample_set[-c(com_ind)]
                          rand_sample_set <- rand_sample_set[1:min_samp]
                          rand_set <- c(rand_set, 
                            rand_sample_set)
                        } else {
                          rand_set <- c(rand_set, 
                            rand_sample_set[1:min_samp])
                          rand_set <- rand_set[order(rand_set)]
                        }
                        rand_set <- na.omit(rand_set)
                        print(rand_set)
                        # EIC.plot(aligned,rows=c(unique(overlapres5ppm$common$index.A)),min.run=runval,min.pres=presval)
                        # apLCMS.EIC.plot(aligned, rows =
                        # c(unique(overlapres5ppm$common$index.A)), colors
                        # = NA, transform = 'none', subset = rand_set,
                        # minrt=NA, maxrt=NA,
                        # min.run=runval,min.pres=presval,
                        # max.spline.time.points = 1000)
                        overlap_res <- overlapres5ppm$common
                        overlap_res <- as.data.frame(overlap_res)
                        dup_mz_ind <- which(duplicated(overlapres5ppm$common$index.A) == 
                          TRUE)
                        if (length(dup_mz_ind) > 0) {
                          overlap_res <- overlap_res[-c(dup_mz_ind), 
                            ]
                        }
                        finalfeatmat <- as.data.frame(finalfeatmat)
                        time.list = finalfeatmat$time[c((overlap_res$index.A))]
                        mz.list = finalfeatmat$mz[c((overlap_res$index.A))]
                        chem.names <- stddata$Name[c((overlap_res$index.B))]
                        
                        eicraw <- getEIC(xset3, groupidx = c(overlap_res$index.A), 
                          rt = "raw")
                        for (i in 1:length(overlap_res$index.A)) {
                          plot(eicraw, xset3, groupidx = i)
                        }
                      }
                      dev.off()
                    }
                    
                  }
                }
            }
        }
    }
    return(aligned_data_list)
}
