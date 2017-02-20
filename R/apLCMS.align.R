apLCMS.align <-
function(cdfloc, apLCMS.outloc,min.run.list=c(3,3), min.pres.list=c(0.3,0.8), minexp=2, mztol=0.00001, alignmztol=NA, alignchrtol=NA,
numnodes=2, run.order.file=NA,subs=NA,filepattern=".cdf",apLCMSmode="untargeted",known_table,match_tol_ppm=5,ref.mz.list=NA,refMZ.mz.diff=10,refMZ.time.diff=NA,target.mz.list = NA)
{
        setwd(cdfloc)
	  cdf.files=list.files(cdfloc,filepattern,ignore.case=TRUE)
	  cdf.files=tolower(cdf.files)
	  dir.create(apLCMS.outloc,showWarnings=FALSE)
	aligned_data_list=new("list")
	pcount=1
	
	
	
        if(is.na(subs[1])==FALSE)
        {
                numsamp=length(subs)
		    cdf.files=cdf.files[subs]
		    
        }
        else
        {
                
                numsamp=length(cdf.files)
        }
        
        if(length(min.run.list)!=length(min.pres.list)){
        	
        	stop("Vectors min.run.list and min.pres.list should be of the same length. eg: min.run.list=c(3,3) and min.pres.list=c(0.3,0.8)")
        }
        for(r in 1:length(min.run.list))
        {
		    runval=min.run.list[r]
		    p=r
                #for(p in 1:length(min.pres.list))
                {
                        features<-new("list")
                        
                        presval=min.pres.list[p]

			par(mfrow=c(2,2))
			fname<-paste("Rplots",runval,presval,".pdf",sep="")
			pdf(fname)
			
			print(minexp)

                        							if(apLCMSmode=="untargeted"){
							aligned<-cdf.to.ftr(cdfloc,subs=subs,min.exp=minexp,min.run=runval,min.pres=presval,mz.tol=mztol, align.mz.tol=alignmztol, align.chr.tol=alignchrtol, n.nodes=numnodes,file.pattern=filepattern)
							}else{
								if(apLCMSmode=="hybrid"){
									
									
aligned<-semi.sup(folder=cdfloc, known.table = known_table, match.tol.ppm = match_tol_ppm, subs=subs,min.exp=minexp,min.run=runval,min.pres=presval,mz.tol=mztol, align.mz.tol=alignmztol, align.chr.tol=alignchrtol, n.nodes=numnodes,file.pattern=filepattern)


									}
								
								}
							
						

					                         
                        
	
			dev.off()
                        finalfeatmat=aligned$final.ftrs
			fname<-paste(apLCMS.outloc, "/apLCMS_aligned", "_pres",  presval, "_run", runval,"_", minexp, "exppostrecovery.Rda", sep="")
			save(aligned,file=fname)

					
			if(is.na(target.mz.list)==FALSE){
			
			eic_fname<-paste(apLCMS.outloc,"/EICrun",runval,"pres",presval,".pdf",sep="")
			
				pdf(eic_fname)
			
				if(is.na(target.mz.list[1,1])==FALSE){
stddata<-target.mz.list

#print(head(stddata))
}else{
	Name<-paste("mz",seq(1,dim(finalfeatmat)[1]),sep="")
	stddata<-cbind(finalfeatmat[,c(1)],Name)
	
	}
				overlapres5ppm<-getVenn(dataA=finalfeatmat, name_a=paste("Expdata",sep=""), dataB=stddata, name_b="Target", mz.thresh = 10, time.thresh=NA,alignment.tool="apLCMS",
xMSanalyzer.outloc=apLCMS.outloc,plotvenn=FALSE)
				if(length(unique(overlapres5ppm$common$index.A))>0)
				{
					setwd(cdfloc)
					num_samples<-dim(finalfeatmat)[2]-4
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
					#EIC.plot(aligned,rows=c(unique(overlapres5ppm$common$index.A)),min.run=runval,min.pres=presval)
					#apLCMS.EIC.plot(aligned, rows = c(unique(overlapres5ppm$common$index.A)), colors = NA, transform = "none",
					#subset = rand_set, minrt=NA, maxrt=NA, min.run=runval,min.pres=presval, max.spline.time.points = 1000)
					#rand_set<-c(1:6)
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
					custom.EIC.plot(aligned, rows = c((overlap_res$index.A)), colors = NA, transform = "none", 
						subset = rand_set, mz.list=mz.list, time.list=time.list,chem.names=chem.names,min.run=runval,min.pres=presval, max.spline.time.points = 1000) 
				}
				dev.off()
			}
                        cnames<-colnames(finalfeatmat[,-c(1:4)])
                        cnames<-tolower(cnames)
                        cnames<-gsub(".cdf", "", cnames)

                        if(is.na(run.order.file)==FALSE)
                        {
                                fileorder=read.table(run.order.file, header=FALSE)
                                fileorder=apply(fileorder,1,tolower)
				cnames=tolower(cnames)
                                ordlist=sapply(1:length(fileorder),function(i){which(cnames==fileorder[i])})
                                ordlist=unlist(ordlist)
                                ordlist=ordlist+4
                                finalfeatmat=finalfeatmat[,c(1:4,ordlist)]
                        }
			fname<-paste(apLCMS.outloc, "/apLCMS_aligned", "_pres", presval, "_run", runval,"_", minexp, "exppostrecovery.txt", sep="")
                        write.table(finalfeatmat,fname,sep="\t",row.names=F)
			
			


			
			aligned_data_list[[pcount]]<-finalfeatmat
			pcount=pcount+1
                }
        }
	return(aligned_data_list)
}
