tic.eval <-
function(dataA,outloc){
	
	dir.create(outloc,showWarnings=FALSE)
	setwd(outloc)
	
	tic<-apply(dataA,2,function(x){
		x<-replace(x,which(x==0),NA)
		return(sum(x,na.rm=TRUE))
	})
	
	
	mean_tic<-mean(tic)

	cv_tic<-100*sd(tic,na.rm=TRUE)/mean(tic,na.rm=TRUE)
	tic_res<-cbind(mean_tic,cv_tic)
	colnames(tic_res)<-c("Average_TIC","CV_TIC")



	main_lab<-paste("Total TIC using all features\n Average TIC=",mean_tic,"\n%CV TIC=",cv_tic,sep="")
	
	
	#tiff("barplot_TIC_using_all_features.tiff",width=2000,height=2000,res=300)
	#pdf("TIC_all_features.pdf")
	barplot(tic,cex.names=0.35,cex.axis=1,main=main_lab,col="orange",cex.main=0.6)
		
	#boxplot(dataA,cex.names=0.35,cex.axis=1,main=main_lab)
	#dev.off()
	
	write.table(tic_res,"TIC_using_all_features.txt",sep="\t",quote=F,col.name=T,row.names=F)
	names(tic)<-c("sample_TIC")
	write.table(tic,"TIC_each_sample_using_all_features.txt",sep="\t",quote=F,col.name=T,row.names=T)
	
	#tiff("boxplot_sampleintensity_using_all_features.tiff",width=2000,height=2000,res=300)
	
}
