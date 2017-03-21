#' getVenn
#' 
#' This utility calls the find.Overlapping.mzs function and generates a Venn
#' diagram showing the extent of overlap between two datasets.
#' 
#' This utility calls the find.Overlapping.mzs and find.Unique.mzs functions
#' and generates a Venn diagram showing the extent of overlap between datasets
#' (up to three).
#' 
#' @param dataA apLCMS or XCMS feature table for dataset A as a data frame.
#' @param name_a Name of dataset A (eg: "SetA").
#' @param dataB apLCMS or XCMS feature table for dataset B as a data frame.
#' @param name_b Name of dataset A (eg: "SetB")
#' @param mz.thresh +/- ppm mass tolerance for m/z matching
#' @param time.thresh Maximum retention time difference (+/-) secs. eg: 30
#' @param alignment.tool "apLCMS" or "XCMS" or "NA". If NA is specified then
#' the first two columns are treated as m/z and retention time, respectively.
#' @param xMSanalyzer.outloc Output folder, eg:
#' "C:/experiment1/xMSanalyzeroutput/"
#' @param use.unique.mz If "TRUE", the function first finds unique features
#' within each set
#' @param plotvenn If "TRUE", the function plots the venn diagram.
#' @return A list is returned.  \item{common }{Row numbers, m/zs, delta
#' retention times of features that are common between the two datasets.}
#' \item{commonA }{Overlapping features in dataset A} \item{uniqueA }{Features
#' that are unique in dataset A} \item{commonB }{Overlapping features in
#' dataset B} \item{uniqueB }{Features that are unique in dataset B}
#' \item{vennCounts}{Output of vennCounts function in limma package}
#' @note Only the unqiue m/zs within each dataset are used to generate Venn
#' diagram.
#' @author Karan Uppal <kuppal2@@emory.edu>
#' @references http://rss.acs.unt.edu/Rdoc/library/limma/html/venn.html
#' @keywords ~Venn
getVenn <- function(dataA, name_a, dataB, name_b, 
    mz.thresh = 10, time.thresh = 30, alignment.tool = NA, 
    xMSanalyzer.outloc, use.unique.mz = FALSE, plotvenn = TRUE) {
    dir.create(xMSanalyzer.outloc, showWarnings = FALSE)
    
    data_a <- as.data.frame(dataA)
    data_b <- as.data.frame(dataB)
    rm(dataA)
    rm(dataB)
    
    ############################################ find.Overlapping.mzs<-function(dataA, dataB,
    ############################################ mz.thresh=10, time.thresh=NA, alignment.tool=NA)
    
    if (use.unique.mz == TRUE) {
        data_a <- find.Unique.mzs.sameset(dataA = data_a, 
            dataB = data_a, mz.thresh = mz.thresh, 
            time.thresh = time.thresh, alignment.tool = alignment.tool)
        data_a <- data_a$uniqueA
        
        # print(dim(data_a))
        
        data_b <- find.Unique.mzs.sameset(dataA = data_b, 
            dataB = data_b, mz.thresh = mz.thresh, 
            time.thresh = time.thresh, alignment.tool = alignment.tool)
        data_b <- data_b$uniqueA
    }
    common <- find.Overlapping.mzs(data_a, data_b, 
        mz.thresh, time.thresh = time.thresh, alignment.tool = alignment.tool)
    
    if (length(common) > 0) {
        commonA <- data_a[c(common$index.A), ]
        commonA <- unique(commonA)
        
        commonB <- data_b[c(common$index.B), ]
        commonB <- unique(commonB)
        
        data_a <- data_a[-c(common$index.A), ]
        data_b <- data_b[-c(common$index.B), ]
        
        
        
        rm_index <- which(data_a$mz %in% common$mz.data.A)
        
        if (length(rm_index) > 0) {
            uniqueA <- data_a[-rm_index, ]
        } else {
            uniqueA <- data_a
        }
        
        rm_index <- which(data_b$mz %in% common$mz.data.B)
        
        if (length(rm_index) > 0) {
            uniqueB <- data_b[-rm_index, ]
        } else {
            uniqueB <- data_b
        }
        
        # uniqueA<-data_a[-c(common$index.A),]
        # uniqueB<-data_b[-c(common$index.B),]
        num_commonA <- length(unique(common$index.A))
        num_commonB <- length(unique(common$index.B))
        
        num_common <- min(num_commonA, num_commonB)[1]
    } else {
        uniqueA <- data_a
        uniqueB <- data_b
        num_common <- 0
        commonA <- {
        }
        commonB <- {
        }
    }
    num_commonA <- num_common
    num_commonB <- num_common
    num_uniqueA <- dim(uniqueA)[1]
    num_uniqueB <- dim(uniqueB)[1]
    
    # print(num_commonA) print(num_uniqueA)
    g1 <- c(seq(1, (num_commonA + num_uniqueA)))
    g2 <- c(seq(1, (num_commonB + num_uniqueB)))
    
    g1[1:num_commonA] = paste("x_", g1[1:num_commonA], 
        sep = "")
    g2[1:num_commonB] = paste("x_", g2[1:num_commonB], 
        sep = "")
    
    if (num_uniqueA > 0) {
        g1[(num_commonA + 1):(num_commonA + num_uniqueA)] = paste("y_", 
            g1[(num_commonA + 1):(num_commonA + num_uniqueA)], 
            sep = "")
    }
    if (num_uniqueA > 0) {
        g2[(num_commonB + 1):(num_commonB + num_uniqueB)] = paste("z_", 
            g2[(num_commonB + 1):(num_commonB + num_uniqueB)], 
            sep = "")
    }
    set1 = as.character(g1)
    set2 = as.character(g2)
    universe <- sort(unique(c(set1, set2)))
    
    Counts <- matrix(0, nrow = length(universe), ncol = 2)
    colnames(Counts) <- c(name_a, name_b)
    for (i in 1:length(universe)) {
        Counts[i, 1] <- universe[i] %in% set1
        Counts[i, 2] <- universe[i] %in% set2
    }
    fname <- paste(xMSanalyzer.outloc, "/Venn", name_a, 
        "_", name_b, "_", mz.thresh, "ppm", time.thresh, 
        "s.pdf", sep = "")
    venn_counts <- vennCounts(Counts)
    if (plotvenn == TRUE) {
        
        pdf(fname)
        vennDiagram(venn_counts)
        dev.off()
    }
    
    return(list(common = common, commonA = commonA, 
        uniqueA = uniqueA, commonB = commonB, uniqueB = uniqueB, 
        vennCounts = venn_counts))
    
}
