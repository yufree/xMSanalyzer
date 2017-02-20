eval.target.mz <-
function(dataA,refMZ,feature.eval.result,mzthresh=10,timethresh=NA,outloc,folderheader=NA,alignment.tool=NA){

feature09<-dataA
if(is.na(refMZ[1,1])==FALSE){
stddata<-refMZ
}else{
	Name<-paste("mz",seq(1,dim(feature09)[1]),sep="")
	stddata<-cbind(feature09[,c(1)],Name)
	
	}
qcresults<-feature.eval.result

rm(dataA)
outloc1<-paste(outloc,"/",folderheader,"targetedeval",mzthresh,"ppm/",sep="")
dir.create(outloc1,showWarnings=FALSE)
setwd(outloc1)

col.names.dataA<-colnames(feature09)

if(is.na(alignment.tool)==FALSE){
	 if(alignment.tool=="apLCMS")
        {
              sample.col.start=5
        }
        else
        {
              if(alignment.tool=="XCMS")
              {
              	
                    sample.col.start=9
                    col.names.dataA[1]="mz"
                    col.names.dataA[2]="time"
		    
                    colnames(data_a)=col.names.dataA
                   
              }
	      
	}}else{
                    #stop(paste("Invalid value for alignment.tool. Please use either \"apLCMS\" or \"XCMS\"", sep=""))
		   
		    col.names.dataA[1]="mz"
		   
		   
		     col.names.dataA[2]="time"
		      sample.col.start=3
		     		      print("Using the 1st columns as \"mz\" and 2nd columsn as \"retention time\"")
		     
		    colnames(feature09)=col.names.dataA
                   
		}
		
par(mfrow=c(2,2))


overlapres5ppm<-getVenn(dataA=feature09, name_a=paste("Expdata",sep=""), dataB=stddata, name_b="Target", mz.thresh = mzthresh, time.thresh=timethresh,alignment.tool="apLCMS",
xMSanalyzer.outloc=outloc1,plotvenn=FALSE)

save(overlapres5ppm,file="overlapRes.Rda")
match5ppmdata<-feature09[overlapres5ppm$common$index.A,]

name_mz<-{}
if(dim(match5ppmdata)[1]<1){
print("No matches found for targeted metabolites.")
return(name_mz);
}


qcresults5ppmdata<-qcresults[overlapres5ppm$common$index.A,]


#m1<-cbind(match5ppmdata[,c(1:(sample.col.start-1))],qcresults5ppmdata[,c(3,7,8)],match5ppmdata[,-c(1:(sample.col.start-1))])

fnames<-paste("../../Stage3b/Targeted_feature_table_",mzthresh,"ppm_filtered.txt",sep="")
#write.table(m1,file=fnames,sep="\t",quote=F,col.name=T,row.names=F)

#write.table(match5ppmdata,file="Matching_data.txt",sep="\t",quote=F,col.name=T,row.names=F)

fnames<-paste("Boxplot_sampleintensity_usingtargetmatches",mzthresh,"ppm.tiff",sep="")


int_data<-log10(match5ppmdata[,-c(1:(sample.col.start-1))]+1)
#int_data<-match5ppmdata[,-c(1:(sample.col.start-1))]
main_lab<-paste("Intensity distribution (log10; all samples) \n of each m/z matching targets (+/- ",mzthresh,"ppm)",sep="")


vennDiagram(circle.col="red",overlapres5ppm$vennCounts,counts.col="blue")

#(fnames,width=2000,height=2000,res=300)
#	boxplot(int_data,cex.names=0.35,cex.axis=1,main=main_lab, ylab="Intensity",xlab="Sample",col="orange",cex.main=0.7)
#	dev.off()

	
	

cv5ppm<-apply(match5ppmdata[,-c(1:(sample.col.start-1))],1,function(x){
	x<-replace(x,which(x==0),NA)
	cvres<-100*sd(x,na.rm=TRUE)/mean(x,na.rm=TRUE)
	return(cvres)
	})
names(cv5ppm)<-round(match5ppmdata[,1],5)

fnames<-paste("TotalCVallsamples_refmz",mzthresh,"ppm.tiff",sep="")

#pdf("Targeted_mz_QC.pdf")
main_name<-paste("Total CV (across all samples) \n of each m/z matching targets (+/- ",mzthresh,"ppm)",sep="")

if(length(cv5ppm)>0){
#tiff(fnames,width=2000,height=2000,res=300)
barplot(cv5ppm,cex.names=0.4,cex.axis=1,main=main_name,col="orange",cex.main=0.7)
}
#dev.off()
fnames<-paste("BarplotTIC_targetmatches",mzthresh,"ppm.tiff",sep="")

main_name<-paste("TIC per sample for matching targets (+/- ",mzthresh,"ppm)",sep="")
#tiff(fnames,width=2000,height=2000,res=300)

tic<-apply(match5ppmdata[,-c(1:(sample.col.start-1))],2,function(x){
	x<-replace(x,which(x==0),NA)
	res<-sum(x,na.rm=TRUE)
	return(res)
})

if(length(tic)>0){
barplot(tic,cex.names=0.4,cex.axis=1,main=main_name,col="orange",cex.main=0.7)
}

fnames<-paste("Pairwiseplot_overall",mzthresh,"ppm.tiff",sep="")
#tiff(fnames,width=2000,height=2000,res=300)
if(length(overlapres5ppm$common$index.A)>0){
plot(cbind(feature09[overlapres5ppm$common$index.A,c(1:2)],qcresults[overlapres5ppm$common$index.A,c(3,7,8)]),main="Pairwise plots of m/z, time, CV, Qscore",cex.main=0.7)
}
#dev.off()



#dev.off()



mean_tic<-mean(tic,na.rm=TRUE)

#cv_tic<-100*sd(tic,na.rm=TRUE)/mean(tic,na.rm=TRUE)

#tic_res<-cbind(mean_tic,cv_tic)

#colnames(tic_res)<-c("Average_TIC","CV_TIC")

#fnames<-paste("TICperfeat_targetmatches",mzthresh,"ppm.txt",sep="")
#write.table(tic_res,file=fnames,sep="\t",quote=F,col.name=T,row.names=F)

tic<-as.data.frame(tic)
names(tic)<-c("sample_TIC")

fnames<-paste("TICpersamp_targetmatches",mzthresh,"ppm.txt",sep="")
write.table(tic,file=fnames,sep="\t",quote=F,col.name=T,row.names=TRUE)


#refMZ<-cbind(refMZ[overlapres5ppm$common$index.B,],match5ppmdata)

#colnames(refMZ)<-c("mz","Name")

name_mz<-cbind(match5ppmdata[,c(1:2)],stddata[overlapres5ppm$common$index.B,],match5ppmdata[,-c(1:2)])

name_mz<-as.data.frame(name_mz)

#name_mz<-merge(refMZ,match5ppmdata,by="mz")

if(folderheader=="raw"){
fnames<-paste("../../Stage3b/",folderheader,"_targetfeaturetable",mzthresh,"ppm_final.txt",sep="")
}else{
	
	if(folderheader=="ComBat"){
	fnames<-paste("../../Stage4b/",folderheader,"_targetfeaturetable",mzthresh,"ppm_final.txt",sep="")
	}
}
#fnames<-paste("Targetmz_matching_data",mzthresh,"ppm.txt",sep="")
write.table(name_mz,file=fnames,sep="\t",quote=F,col.name=T,row.names=F)


#outloc1<-paste(outloc1,"target_barplots/",sep="")
#dir.create(outloc1,showWarnings=FALSE)
#setwd(outloc1)

match5ppmdata<-name_mz
match5ppmdata<-as.data.frame(match5ppmdata)
match5ppmdata_int<-match5ppmdata[,-c(1:(2+dim(refMZ)[2]))]
match5ppmdata_int<-apply(match5ppmdata_int,2,as.numeric)

print(head(match5ppmdata))
print(head(match5ppmdata_int))
print(dim(match5ppmdata))
print(dim(match5ppmdata_int))
if(length(overlapres5ppm$common$index.B)>0){
for(i in 1:dim(match5ppmdata_int)[1]){
	fname<-paste("mz",round(match5ppmdata[i,1],5),".tiff",sep="")
#	tiff(fname,width=2000,height=2000,res=300)
delta_ppm<-10^6*(abs(match5ppmdata[i,1]-match5ppmdata[i,3]))/match5ppmdata[i,3]
delta_ppm<-round(delta_ppm,2)

#int_val<-as.numeric(int_val)
if(length(overlapres5ppm$common$index.B)>1){
cur_vec<-t(match5ppmdata_int[i,])
}else{
cur_vec<-t(match5ppmdata_int[i])
}
	plot(x=t(cur_vec),main=paste("mz: ",round(match5ppmdata[i,1],5),"; time:",round(match5ppmdata[i,2],0), "(s) ;\n",match5ppmdata$Name[i], " (delta=",delta_ppm," ppm)",sep=""),ylab="Intensity",xlab="Sample index",cex=0.6,cex.main=0.7,col="brown")
	#dev.off()
}

}
#dev.off()

return(name_mz)
}
