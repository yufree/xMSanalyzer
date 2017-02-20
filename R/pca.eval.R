pca.eval <-
function(X,samplelabels,filename=NA,ncomp=5,center=TRUE,scale=TRUE,legendlocation="topright",legendcex=0.5,outloc=getwd()){
	
	
	X<-as.matrix(X)
	print("Performing PCA")
	batchlabels<-samplelabels
	metabpcaresultlog2allmetabs5pcs<-pca(X,ncomp=ncomp,center=center,scale=scale)
    
    #metabpcaresultnotransform10pcsallmetabs<-metabpcaresult
    result<-metabpcaresultlog2allmetabs5pcs
    
    r1<-100*result$sdev/sum(result$sdev)
    r1<-round(r1,2)
    
    col_vec<-c("mediumpurple4","mediumpurple1","blueviolet","darkblue","blue","cornflowerblue","cyan4","skyblue",
    "darkgreen", "seagreen1", "green","yellow","orange","pink", "coral1", "palevioletred2",
    "red","saddlebrown","brown","brown3","white","darkgray","aliceblue",
    "aquamarine","aquamarine3","bisque","burlywood1","lavender","khaki3","black")
    
    #col_vec<-col_vec[sample(1:length(col_vec),length(col_vec))]
    
    t1<-table(samplelabels)
    l1<-levels(samplelabels)
    col_all=topo.colors(256)
    
    dir.create(outloc,showWarnings=FALSE)
    setwd(outloc)
    print(paste("Generating PCA plots to evaluate batch-effect in ",filename," data",sep=""))
    
    fname<-paste("PCA_batcheffecteval",filename,".tiff",sep="")
    ## 1) raw data
    #tiff(fname,res=300, width=2000,height=2000)
    col <- rep(col_vec[1:length(t1)], t1)
    
    #col<-rep(col_all[1:length(l1)],t1)
    ## Choose different size of points
    cex <- rep(1, dim(X)[1])
    
    ## Choose the form of the points (square, circle, triangle and diamond-shaped
    #pch <- replace(pch,which(pch==2),21)
    #pch <- replace(pch,which(pch==3),17)
    
    pch <- rep(15,dim(X)[1])
    
    ## comp is the two dimensions you want to display
    ## ind.names is whether you want labels on each points
    ## rep.space determines the subspace to project the individuals ("X-variate",
    ## "Y-variate" or "XY-variate")
    #plotIndiv(result, comp = c(1,2), ind.names = FALSE, rep.space = "X-variate", col = col, cex = cex, pch = pch, X.label="PC1",Y.label="PC2")
    
    pca_res<-suppressWarnings(try(plotIndiv(result, comp = c(1,2), ind.names = FALSE, col = col, cex = cex, pch = pch, X.label=paste("PC1 (",r1[1],"% variation)",sep=""),Y.label=paste("PC2 (",r1[2],"% variation)",sep="")),silent=TRUE))
   
	if(is(pca_res, "try-error")){
		
		print("PCA plot could not be generated.")
	}else{ 
    col<-col_vec[1:length(t1)]
    cex <- rep(0.4,length(t1))
    pch <- rep(15,length(t1))
    
    ## The first two parameters are the x and y coordinates of the legend on the graph
    ## The third one is the text of the legend
    legend(legendlocation, l1, col = col,pch = pch, pt.cex = cex, title = "Batch #", cex=legendcex)
    
   # dev.off()
    }

}
