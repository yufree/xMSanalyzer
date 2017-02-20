getVennmultiple <-
function(dataA,name_a, dataB,name_b,dataC,name_c,mz.thresh=10,time.thresh=30,alignment.tool=NA, xMSanalyzer.outloc, use.unique.mz=FALSE,plotvenn=TRUE)
{
	dir.create(xMSanalyzer.outloc,showWarnings=FALSE)
	
	 data_a<-as.data.frame(dataA)
        data_b<-as.data.frame(dataB)
	 data_c<-as.data.frame(dataC)
	
	rm(dataA)
	rm(dataB)
	rm(dataC)
	
	data_a<-as.data.frame(data_a)
	data_b<-as.data.frame(data_b)
	data_c<-as.data.frame(data_c)

	if(use.unique.mz==TRUE){
	    data_a<-find.Unique.mzs.sameset(data_a,data_a,mz.thresh=mz.thresh,time.thresh=time.thresh,alignment.tool=alignment.tool)
	    data_a<-data_a$uniqueA
	    
	data_b<-find.Unique.mzs.sameset(data_b,data_b,mz.thresh=mz.thresh,time.thresh=time.thresh,alignment.tool=alignment.tool)
	 data_b<-data_b$uniqueA
	    
	    data_c<-find.Unique.mzs.sameset(data_c,data_c,mz.thresh=mz.thresh,time.thresh=time.thresh,alignment.tool=alignment.tool)
	     data_c<-data_c$uniqueA
	}
	commonAB<-find.Overlapping.mzs(data_a,data_b,mz.thresh,time.thresh,alignment.tool=alignment.tool)
	commonAC<-find.Overlapping.mzs(data_a,data_c,mz.thresh,time.thresh,alignment.tool=alignment.tool)
	commonBC<-find.Overlapping.mzs(data_b,data_c,mz.thresh,time.thresh,alignment.tool=alignment.tool)
	
	#if(length(commonAB)>0)
	{
	data_ab<-data_a[commonAB$index.A,]
	
	data_ab<-as.data.frame(data_ab)
	
	commonABC<-find.Overlapping.mzs(data_ab,data_c,mz.thresh,time.thresh,alignment.tool=alignment.tool)
	}
	data_ba<-data_b[commonAB$index.B,]
	data_ba<-as.data.frame(data_ba)
	
	commonBAC<-find.Overlapping.mzs(data_ba,data_c,mz.thresh,time.thresh=NA,alignment.tool=alignment.tool)
	
	#get unique A
	rm_index<-which(data_a$mz%in%commonAB$mz.data.A)
	
	rm_index<-c(rm_index,which(data_a$mz%in%commonAC$mz.data.A))
	
	if(length(rm_index)>0){
	uniqueA<-data_a[-rm_index,]
	uniqueA<-as.data.frame(uniqueA)
	}
	

	rm_index<-which(data_b$mz%in%commonAB$mz.data.B)
	
	rm_index<-c(rm_index,which(data_b$mz%in%commonBC$mz.data.A))
	
	#rm_index<-which(uniqueB$mz%in%commonBC$mz.data.A)

	if(length(rm_index)>0){
	uniqueB<-data_b[-rm_index,]
	uniqueB<-as.data.frame(uniqueB)
	}	
	
	
	rm_index<-which(data_c$mz%in%commonAC$mz.data.B)

	rm_index<-c(rm_index,which(data_c$mz%in%commonBC$mz.data.B))

	if(length(rm_index)>0){
	uniqueC<-data_c[-rm_index,]
	uniqueC<-as.data.frame(uniqueC)
	}
	
	
	
	num_commonAB<-length(unique(commonAB$mz.data.A))-length(which(unique(commonAB$mz.data.A)%in%unique(commonABC$mz.data.A)))
	
	num_commonBC<-length(unique(commonBC$mz.data.A))-length(which(unique(commonBC$mz.data.A)%in%unique(commonBAC$mz.data.A)))
	
	num_commonAC<-length(unique(commonAC$mz.data.A))-length(which(unique(commonAC$mz.data.A)%in%unique(commonABC$mz.data.A)))
	
	num_commonCB<-length(unique(commonBC$mz.data.B))-length(which(unique(commonBC$mz.data.B)%in%unique(commonBAC$mz.data.B)))
	
	num_commonCA<-length(unique(commonAC$mz.data.B))-length(which(unique(commonAC$mz.data.B)%in%unique(commonABC$mz.data.B)))
	
	
	num_commonABC<-min(length(unique(commonABC$mz.data.A)),length(unique(commonABC$mz.data.B)))
	
	num_uniqueA<-dim(uniqueA)[1]
	num_uniqueB<-dim(uniqueB)[1]
	num_uniqueC<-dim(uniqueC)[1]
	
	g1 <-paste("a",seq(num_commonAB+num_commonAC+num_commonABC+num_uniqueA),sep="")

	g2<-paste("b",seq(num_commonAB+num_commonBC+num_commonABC+num_uniqueB),sep="")
	
	g3<-paste("c",seq(num_commonCA+num_commonCB+num_commonABC+num_uniqueC),sep="")

	#x: AB; w:AC; v:BC;u:ABC
	if(num_commonAB>0)
	{
		g1[1:num_commonAB]=paste("x_",seq(1,num_commonAB),sep="")
		
		g2[1:num_commonAB]=paste("x_",seq(1,num_commonAB),sep="")
	}
	
	if(num_commonAC>0)
	{
		g1[(num_commonAB+1):(num_commonAB+num_commonAC)]=paste("w_",seq(1,num_commonAC),sep="")
		
		g3[1:num_commonAC]=paste("w_",seq(1,num_commonAC),sep="")
		
	}
	
	if(num_commonBC>0)
	{
		g2[(num_commonAB+1):(num_commonAB+num_commonBC)]=paste("v_",seq(1,num_commonBC),sep="")
	
		g3[(num_commonAC+1):(num_commonAC+num_commonBC)]=paste("v_",seq(1,num_commonBC),sep="")
	}
	
	g1[(num_commonAB+num_commonAC+1):(num_commonAB+num_commonAC+num_commonABC)]=paste("u_",seq(1,num_commonABC),sep="")
	g2[(num_commonAB+num_commonBC+1):(num_commonAB+num_commonBC+num_commonABC)]=paste("u_",seq(1,num_commonABC),sep="")
	g3[(num_commonAC+num_commonBC+1):(num_commonAC+num_commonBC+num_commonABC)]=paste("u_",seq(1,num_commonABC),sep="")
	

	set1=as.character(g1)
	set2=as.character(g2)
	set3=as.character(g3)
	universe <- sort(unique( c(set1,set2,set3)))
	Counts <- matrix(0, nrow=length(universe), ncol=3)
	colnames(Counts) <- c(name_a, name_b,name_c)
	for (i in 1:length(universe))
	{
		Counts[i,1] <- universe[i] %in% set1
		Counts[i,2] <- universe[i] %in% set2
		Counts[i,3] <- universe[i] %in% set3
	}
	venn_counts<-vennCounts(Counts)
	fname<-paste(xMSanalyzer.outloc,"/Venn", name_a,"_",name_b,"_",name_c,"_",mz.thresh,"ppm",time.thresh,"s.pdf",sep="")
	
	if(plotvenn==TRUE){
		pdf(fname)
		vennDiagram(venn_counts)
		dev.off()
	}
	return(list("commonABC"=commonABC,"uniqueA"=uniqueA,"uniqueB"=uniqueB,"uniqueC"=uniqueC,"commonAB"=commonAB, 
	"commonBC"=commonBC,"commonAC"=commonAC,"vennCounts"=venn_counts))

}
