check.mz.in.replicates <-
function(dataA,min.samps=2,min.reps=2,num_replicates=3)
{
        mean_replicate_difference<-{}
        sd_range_duplicate_pairs<-{}
	
	curdata<-dataA
	rm(dataA)

	#curdata<-curdata[1:10,]
        numfeats=dim(curdata)[1]
        numsamp=dim(curdata)[2]
        rnames<-colnames(curdata)
        rnames<-gsub(".cdf", "", rnames, ignore.case=TRUE)
	quantcolnames=c("min", "first_quartile", "median", "mean", "third_quartile", "max")
	

        
                newrow={}
                finalmat={}
		
		
		cl<-makeSOCKcluster(10)
			
		
		
		clusterExport(cl, "check.mz.in.replicates.child") 
		
		cv.res<-parApply(cl,curdata,1,check.mz.in.replicates.child,min.samps=min.samps,min.reps=min.reps,num_replicates=num_replicates)
		print("done")
	dim(cv.res)=dim(matrix(nrow=1,ncol=numfeats))	
	print("done")
	stopCluster(cl)
        #final_set<-as.data.frame(cv.res)
	
        #rownames(final_set)=NULL
	#final_set<-apply(final_set,2,as.numeric)
	#final_set<-as.data.frame(t(final_set))
	#colnames(final_set)<-quantcolnames
        final_set<-curdata[which(cv.res==1),]
	#final_set<-cbind(curdata_mz_rt_info,final_set)
	return(final_set)
}
