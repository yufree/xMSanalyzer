merge.Results <-
function(dataA, dataB, feat.eval.A,feat.eval.B,max.mz.diff=15, max.rt.diff=300,merge.eval.pvalue=0.05,alignment.tool="apLCMS", numnodes=1, mult.test.cor=FALSE,
mergecorthresh=0.7,missingvalue=0)
{
	dataA<-unique(dataA)
	dataB<-unique(dataB)
	

	  if(alignment.tool=="apLCMS")
        {
              sample.col.start=5
	      
	      dataA<-cbind(dataA,feat.eval.A[,c(sample.col.start:(dim(feat.eval.A)[2]))])
	      dataB<-cbind(dataB,feat.eval.B[,c(sample.col.start:(dim(feat.eval.B)[2]))])
	curdata<-rbind(dataA,dataB)
	curdata<-curdata[order(curdata$mz),]
	curdata<-as.data.frame(curdata)
	
			
	      npeaks<-apply(curdata[,sample.col.start:(dim(curdata)[2]-8)],1,countpeaks)
	      curdatatemp<-cbind(curdata[,c(1:4)],npeaks)
	      curdata<-cbind(curdatatemp,curdata[,sample.col.start:dim(curdata)[2]])
	      sample.col.start=6
	     
	      curdata<-as.data.frame(curdata)
	     
        }
        else
        {
              if(alignment.tool=="XCMS")
              {
                    sample.col.start=9
		    
				    cnames<-colnames(dataA)
                    cnames[1]="mz"

                    cnames[4]="time"
                    colnames(dataA)=cnames
                    colnames(dataB)=cnames
		    
		       	dataA<-cbind(dataA,feat.eval.A[,c(sample.col.start:(dim(feat.eval.A)[2]))])
				dataB<-cbind(dataB,feat.eval.B[,c(sample.col.start:(dim(feat.eval.B)[2]))])
				curdata<-rbind(dataA,dataB)
			
				curdata<-curdata[order(curdata$mz),]
				curdata<-as.data.frame(curdata)
			#curdata<-curdata[-which(curdata$min==curdata$max),]
              }
              else
              {
                    stop(paste("Invalid value for alignment.tool. Please use either \"apLCMS\" or \"XCMS\"", sep=""))
              }
        }
	
        #dataA<-dataA[order(dataA$mz),]
	    #dataB<-dataB[order(dataB$mz),]
	
	
	sub_data_a<-new("list")
	sub_data_b<-new("list")
	lindex<-1
	
	min_mz<-min(curdata[,1])
	max_mz<-max(curdata[,1])

	mz_group=10

	mzdefect<-1*((curdata$mz-floor(curdata$mz)))

	#curdata<-cbind(mzdefect,curdata)
	
	curdata<-as.data.frame(curdata)
	curdata<-unique(curdata)
    max_rt_thresh<-max(curdata$time,na.rm=TRUE)
	
	#d1<-density(mzdefect,bw="nrd",from=min(mzdefect),to=(0.01+max(mzdefect)))

	#sub_data_a<-sapply(list(myData1=curdata),function(x)  split(x,cut(curdata$mz,breaks=seq(min(curdata$mz)-0.001,max(curdata$mz)+d1$bw,d1$bw))))
	#d1$bw=d1$bw*2
	
	#sub_data_a<-sapply(list(myData1=curdata),function(x)  split(x,cut(curdata$mzdefect,breaks=seq(0,1,0.03))))
 
# i<-2
 #j<-3
 #10^6*abs((850+max(sub_data_a[[i]][,1]))-(850+min(sub_data_a[[j]][,1])))/(850+min(sub_data_a[[j]][,1]))
	
	
	
	
	#lindex<-length(sub_data_a)
	diff_mz_num=1
	 #Step 1 Group features by m/z
        mz_groups<-lapply(1:dim(curdata)[1],function(j){
                                commat={}
                                diffmz=new("list")
                                ppmb=(500)*(curdata$mz[j]/1000000)
                                getbind_same<-which(abs(curdata$mz-curdata$mz[j])<=ppmb)
                                diffmz[[diff_mz_num]]=getbind_same #dataA[getbind_same,]
                                diff_mz_num=diff_mz_num+1
                                return(diffmz)
          })
	  	  
	  	  options(warn=-1)
	  	  
	  del_list<-{}
	  for(m in 1:length(mz_groups)){
	  
	   if((m%in%del_list)==FALSE){
		for(n in (m+1):length(mz_groups)){
		
			if(n>length(mz_groups)){
				break;
			}
			com1<-intersect(mz_groups[[m]][[1]],mz_groups[[n]][[1]])
			if(length(com1)>0){
			
				mz_groups[[m]][[1]]<-c(mz_groups[[m]][[1]],mz_groups[[n]][[1]])
				del_list<-c(del_list,n)
			}
		}
		
		mz_groups[[m]][[1]]<-unique(mz_groups[[m]][[1]])
		
		}
	  }
	  if(length(del_list)>0){
		 mz_groups<-mz_groups[-del_list]
		}
	 
	  
	mz_groups<-unique(mz_groups)
	

	sub_data_a<-lapply(1:length(mz_groups),function(j){
                                commat={}
                               # diffmz=new("list")
                                getbind_same=mz_groups[[j]][[1]]
                                diffmz=curdata[getbind_same,]
                                
                                return(diffmz)
          })
          rm(mz_groups)
	lindex<-length(sub_data_a)
	
	sub_data_a<-unique(sub_data_a)
	
	
	
	if(lindex>2){
		cl<-makeSOCKcluster(numnodes)
	
		
			clusterEvalQ(cl, "merge.Results.child.ttest")
			clusterEvalQ(cl, "compare_intensities_ttest")
			
			merge.res<-parLapply(cl,sub_data_a,merge.Results.child.ttest,max.mz.diff=max.mz.diff, max.rt.diff=max.rt.diff, merge.eval.pvalue=merge.eval.pvalue,alignment.tool=alignment.tool, mult.test.cor=mult.test.cor,mergecorthresh=mergecorthresh,col.rm.index=NA,missingvalue=missingvalue)
		
		if(FALSE){
		merge.res<-new("list")
		
		for(lsub in 1:length(sub_data_a)){
			
			tempres<-merge.Results.child.ttest(dataA=sub_data_a[[lsub]],max.mz.diff=max.mz.diff, max.rt.diff=max.rt.diff, merge.eval.pvalue=merge.eval.pvalue,alignment.tool=alignment.tool, mult.test.cor=mult.test.cor,mergecorthresh=mergecorthresh,col.rm.index=NA,missingvalue=missingvalue)
			
			merge.res[[lsub]]<-tempres
			
		}
		}	
	
		stopCluster(cl)
		
		
        #save(merge.res,file="mergeres.Rda")
		
		
	end_ind<-dim(curdata)[2]-8
	l1<-lapply(1:length(merge.res),function(x){nrow(merge.res[[x]])})
	l2<-unlist(l1)
	
	if(length(which(l2<1))>0){
	merge.res<-merge.res[-which(l2<1)]	
	}
	
	l1<-lapply(1:length(merge.res),function(x){is.na(merge.res[[x]][1])[1]})
	l2<-unlist(l1)
	
	if(length(which(l2==TRUE))>0){
	merge.res<-merge.res[-which(l2==TRUE)]	
	}
	
	final.res={}
	
	
	
	#final.res<-ldply(merge.res,rbind)
	#final.res={}
	for(s in 1:length(merge.res))
	{	
		tempd<-as.data.frame(merge.res[[s]])
		final.res<-rbind(final.res,tempd)
		
	}
	
	
	
	final.res<-unique(final.res)
	
	}else{
		
final.res<-merge.Results.child.ttest(curdata,max.mz.diff=max.mz.diff, max.rt.diff=max.rt.diff, merge.eval.pvalue=merge.eval.pvalue,alignment.tool=alignment.tool, mult.test.cor=mult.test.cor,mergecorthresh=mergecorthresh,col.rm.index=NA,missingvalue=missingvalue)
		}
	
	final.res<-as.data.frame(final.res)
	#final.res<-na.omit(final.res)
	final.res<-unique(final.res)
	
	
	final.res<-final.res[order(final.res$mz),]
	
	
	
	#if(FALSE)
	{
	curdata<-final.res
	curdata<-as.data.frame(curdata)
	
	
    #system.time(global_cor<-WGCNA::cor(t(curdata[,c(sample.col.start:end_ind)]),nThreads=numnodes,method="spearman",use = 'p'))

	#count_cor<-lapply(1:dim(curdata)[1],function(j){length(which(global_cor[j,]>=0.9))})
	
	#count_cor<-lapply(1:dim(curdata)[1],function(j){diff_rt<-abs(curdata[j,2]-curdata[,2]); length(which(global_cor[j,]>=0.9 & diff_rt<10))})
	
	
	count_cor<-lapply(1:dim(curdata)[1],function(j){diff_mz<-10^6*abs(curdata[j,1]-curdata[,1])/curdata[j,1]; diff_rt<-abs(curdata[j,2]-curdata[,2]); length(which(diff_mz<max.mz.diff & diff_rt<90))})
	
	count_cor<-unlist(count_cor)
	
	#print("length of count_cor")
	#print(length(count_cor))
	
	
	if(length(which(count_cor<=1))>0){
	non_cor_data<-curdata[which(count_cor<=1),]
	curdata<-curdata[-which(count_cor<=1),]
	#global_cor<-global_cor[-which(count_cor<=1),-which(count_cor<=1)]
	

	
	
	}

	if(dim(curdata)[1]>0){
	#print("dim of curdata 2nd round check")
	#print(dim(curdata))
	final.res2<-merge.Results.child.ttest(curdata,max.mz.diff=max.mz.diff, max.rt.diff=max.rt.diff, merge.eval.pvalue=merge.eval.pvalue,alignment.tool=alignment.tool, mult.test.cor=mult.test.cor,mergecorthresh=mergecorthresh,col.rm.index=NA,missingvalue=missingvalue)
	#final.res2<-na.omit(final.res2)
	final.res2<-unique(final.res2)
    
	final.res<-rbind(final.res2,non_cor_data)
	final.res<-as.data.frame(final.res)
	final.res<-final.res[order(final.res$mz),]
	

	}
    
    curdata<-final.res
    #rm(final.res)
    
    curdata<-replace(as.matrix(curdata),which(curdata==0),NA)
    
    system.time(global_cor<-WGCNA::cor(t(curdata[,c(sample.col.start:end_ind)]),nThreads=numnodes,use = 'p'))
    
    
    count_cor<-lapply(1:dim(curdata)[1],function(j){diff_mz<-10^6*abs(curdata[j,1]-curdata[,1])/curdata[j,1]; diff_rt<-abs(curdata[j,2]-curdata[,2]); length(which(diff_mz<max.mz.diff & global_cor[j,]>=0.9))})
    
    count_cor<-unlist(count_cor)
    
    #print("length of count_cor")
    #print(length(count_cor))
    
    
    if(length(which(count_cor<=1))>0)
    {
        non_cor_data<-curdata[which(count_cor<=1),]
        curdata<-curdata[-which(count_cor<=1),]
        #global_cor<-global_cor[-which(count_cor<=1),-which(count_cor<=1)]
        
        
    }
    
    if(dim(curdata)[1]>0){
        #print("dim of curdata 2nd round check")
        #print(dim(curdata))
        
        curdata<-as.data.frame(curdata)
       
        
        final.res2<-merge.Results.child.ttest(curdata,max.mz.diff=max.mz.diff, max.rt.diff=max_rt_thresh, merge.eval.pvalue=merge.eval.pvalue,alignment.tool=alignment.tool, mult.test.cor=mult.test.cor,mergecorthresh=mergecorthresh,col.rm.index=NA,missingvalue=NA)
        
        rm(curdata)
        #final.res2<-na.omit(final.res2)
        final.res2<-unique(final.res2)
        
        print("here2")
        final.res2<-replace(as.matrix(final.res2),which(is.na(final.res2)==TRUE),0)
        
        final.res<-rbind(final.res2,non_cor_data)
        final.res<-as.data.frame(final.res)
        final.res<-final.res[order(final.res$mz),]
        
        
    }

    
    
    
    
    
    
    
	
	}
	
	if(is.na(missingvalue)==FALSE){
		
			final.res<-replace(as.matrix(final.res),which(is.na(final.res)==TRUE),missingvalue)
			
		}

	final.res<-as.data.frame(final.res)
	
	print("dim of final res")
	print(dim(final.res))
	
	options(warn=0)
	return(final.res)
	#return(sub_data_a)
}
