XCMS.align.centWave <-
function(cdfloc, XCMS.outloc,ppm.list=c(10), mz.diff.list=c(-0.001), sn.thresh.list=c(10), prefilter.list=c(3,100), bw.val=c(10),groupval.method="medret", 
step.list=c(0.1),max=50,minfrac.val=0.5, minsamp.val=1, mzwid.val=0.25, sleep.val=0, run.order.file=NA,subs=NA, retcor.method="obiwarp",retcor.family="symmetric", retcor.plottype="deviation", peakwidth=c(20,50),nSlaves=2,target.mz.list = NA)
{

        setwd(cdfloc)
	dir.create(XCMS.outloc,showWarnings=FALSE)
        cdf_files=list.files(cdfloc,".cdf|.mzxml|mXML",ignore.case=TRUE)
	aligned_data_list=new("list")
	pcount=1
	
	
	
        if(is.na(subs[1])==FALSE)
        {
                cdf_files=cdf_files[subs]
                numsamp=length(subs)


        }
        else
        {

                numsamp=length(cdf_files)
        }

        for(t in sn.thresh.list)
        {
                for(p in ppm.list)
                {
                        for(m in mz.diff.list)
                        {
                                #for(maxl in max.list)
                                for(s in step.list)
				{
					
					xset=xcmsSet(cdf_files, method="centWave", ppm=p, snthresh=t,mzdiff=m,peakwidth=peakwidth,prefilter=prefilter.list,
					integrate=1, verbose.columns=TRUE,fitgauss=FALSE, nSlaves=nSlaves)
					
					####group.density
                                        xset<-group(xset)
					
					if(retcor.method=="loess"){
					
					#####retention time correction using loess; mdevden
                                        xset2 <- retcor(xset, method=retcor.method,family = retcor.family, plottype = retcor.plottype)
					}else{
						if(retcor.method=="obiwarp"){
						###retention time correction using obiwarp method
						xset2 <- retcor(xset, method=retcor.method,plottype = retcor.plottype, profStep=s)
						}
						else{
						stop("Please enter obiwarp or loess as retention time correction method.")
						}
					}
					#retcor(object, method="obiwarp", plottype = c("none", "deviation"), profStep=1, center=NULL, col = NULL, 
					#ty = NULL, response=1, distFunc="cor_opt", gapInit=NULL, gapExtend=NULL, factorDiag=2, factorGap=1, localAlignment=0, initPenalty=0) 
					#findPeaks.centWave(object, ppm=25, peakwidth=c(20,50), snthresh=10, prefilter=c(3,100), 
					#mzCenterFun="wMean", integrate=1, mzdiff=-0.001, fitgauss=FALSE, scanrange= numeric(), noise=0, sleep=0, verbose.columns=FALSE, ROI.list=list()) 
					#retcor(object, method="obiwarp", plottype = c("none", "deviation"), profStep=1, center=NULL, col = NULL, 

					#ty = NULL, response=1, distFunc="cor_opt", gapInit=NULL, gapExtend=NULL, factorDiag=2, factorGap=1, localAlignment=0, initPenalty=0) 
					
					
					
					for(b in bw.val){
                                        ###  Group peaks together across samples, set bandwitdh, change important m/z parameters here
                                        ###  Syntax: group(object, bw = 30, minfrac = 0.5, minsamp= 1,  mzwid = 0.25, max = 5, sleep = 0)
                                       
                                        xset2 <- group(xset2, bw=b, minfrac=minfrac.val, minsamp=minsamp.val,  mzwid=mzwid.val, max=max, sleep=sleep.val)
                                        #xset3=fillPeaks(xset2)
                                       print(xset2)
       					xset3=suppressWarnings(fillPeaks(xset2,method="chrom"))                                 
					print(xset3)
					 
					print("Getting sample intensities")
                                        finalfeatmat={}
					
					if (dim(xset3@groups)[1] > 0) {
					     groupmat <- groups(xset3)
					     #group_intmat<-groupval(xset3, method = groupval.method, intensity = "into") 

					     group_intmat<-groupval(xset3,groupval.method,"into")
					     finalfeatmat<-cbind(groupmat,group_intmat)
					     finalfeatmat<-as.data.frame(finalfeatmat,row.names = NULL)
					     
					  } else{
						if (length(xset3@sampnames) == 1)
						{
							finalfeatmat  <- xset3@peaks
						}else 
						{
							stop ('First argument must be a xcmsSet with group information or contain only one sample.')
						 }
					}

					fname=paste(XCMS.outloc, "/XCMS_centwave","_snthresh", t,"_step", s, "_mzdiff", m,"_max",max,"_bw",b,"_ppm",p,".Rda", sep="")
					save(xset3,file=fname)

                                        fname=paste(XCMS.outloc, "/XCMS_centwave","_snthresh", t,"_step", s, "_mzdiff", m,"_max",max,"_bw",b,"_ppm",p,".txt", sep="")
                                        cnames=colnames(finalfeatmat)
                                        cnames[1]="mz"
                                        cnames[4]="time"
                                        colnames(finalfeatmat)=c(cnames[1:8],cdf_files)
                                        cnames<-tolower(cnames)
                                        cnames<-gsub(".cdf", "", cnames)

                                        if(is.na(run.order.file)==FALSE)
                                        {
                                                fileorder=read.table(run.order.file, header=FALSE)
                                                fileorder=apply(fileorder,1,tolower)
                                                ordlist=sapply(1:length(fileorder),function(i){which(cnames==fileorder[i])})
                                                ordlist=unlist(ordlist)
                                                ordlist=ordlist+8
                                                finalfeatmat=finalfeatmat[,c(1:8,ordlist)]
                                        }

                                        write.table(finalfeatmat,file=fname,sep="\t",row.names=FALSE)
					aligned_data_list[[pcount]]<-finalfeatmat
					pcount=pcount+1
					if(is.na(target.mz.list)==FALSE){

        fname=paste(XCMS.outloc, "/EICmatchedFilter","_thresh", t,"_step", s, "_mzdiff", m,"_max",max,"_bw",b,".pdf", sep="")
                                        pdf(fname)
                                if(is.na(target.mz.list[1,1])==FALSE){
stddata<-target.mz.list

#print(head(stddata))
}else{
        Name<-paste("mz",seq(1,dim(finalfeatmat)[1]),sep="")
        stddata<-cbind(finalfeatmat[,c(1)],Name)

        }
overlapres5ppm<-getVenn(dataA=finalfeatmat, name_a=paste("Expdata",sep=""), dataB=stddata, name_b="Target", mz.thresh = 10, time.thresh=NA,alignment.tool="XCMS",
xMSanalyzer.outloc=XCMS.outloc,plotvenn=FALSE)
                                if(length(unique(overlapres5ppm$common$index.A))>0)
                                {
                                        setwd(cdfloc)
                                        num_samples<-dim(finalfeatmat)[2]-8
                                        min_samp<-6
                                        if(num_samples<min_samp){
                                                min_samp<-num_samples
                                        }
                                        rand_sample_set<-sample(size=num_samples,x=num_samples,replace=FALSE)
                                        rand_set<-c(1:min_samp,(num_samples-2):num_samples)
                                        com_ind<-which(rand_sample_set%in%rand_set)
                                        if(length(com_ind)>0){
                                        rand_sample_set<-rand_sample_set[-c(com_ind)]
                                        rand_sample_set<-rand_sample_set[1:min_samp]
                                        rand_set<-c(rand_set,rand_sample_set)
                                        }else{
                                                rand_set<-c(rand_set,rand_sample_set[1:min_samp])
                                                rand_set<-rand_set[order(rand_set)]
                                        }
                                        rand_set<-na.omit(rand_set)
                                        print(rand_set)
                                        overlap_res<-overlapres5ppm$common
                                        overlap_res<-as.data.frame(overlap_res)
                                        dup_mz_ind<-which(duplicated(overlapres5ppm$common$index.A)==TRUE)
                                        if(length(dup_mz_ind)>0){
                                                overlap_res<-overlap_res[-c(dup_mz_ind),]
                                        }
                                        finalfeatmat<-as.data.frame(finalfeatmat)
					time.list=finalfeatmat$time[c((overlap_res$index.A))]
                                        mz.list=finalfeatmat$mz[c((overlap_res$index.A))]
                                        chem.names<-stddata$Name[c((overlap_res$index.B))]

                                        eicraw <- getEIC(xset3, groupidx = c(overlap_res$index.A), rt = "raw")
                                        for(i in 1:length(overlap_res$index.A)){
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
