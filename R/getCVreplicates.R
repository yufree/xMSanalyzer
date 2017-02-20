getCVreplicates <-
function(curdata,alignment.tool,numreplicates, min.samp.percent=0.6,impute.bool=TRUE,missingvalue=0)
{
        mean_replicate_difference<-{}
        sd_range_duplicate_pairs<-{}
	
        if(alignment.tool=="apLCMS")
        {
              col_end=4
        }
        else
        {
              if(alignment.tool=="XCMS")
              {
                    col_end=8
              }
              else
              {
			  col_end=2
		     print("**Using the first two columns as mz and retention time for PID calculation**")
                    #stop(paste("Invalid value for alignment.tool. Please use either \"apLCMS\" or \"XCMS\"", sep=""))
              }
        }
        
        curdata_mz_rt_info=curdata[,c(1:col_end)]
        curdata=curdata[,-c(1:col_end)]
	
	#curdata<-curdata[1:10,]
        numfeats=dim(curdata)[1]
        numsamp=dim(curdata)[2]
        rnames<-colnames(curdata)
        rnames<-gsub(".cdf", "", rnames, ignore.case=TRUE)
	quantcolnames=c("min", "first_quartile", "median", "mean", "third_quartile", "max","sampleCount")
	if(impute.bool==TRUE)
	{
		min.samp.percent=min.samp.percent
	}
	else
	{
		min.samp.percent=1
	}

        
                newrow={}
                finalmat={}
		
		
		cl<-makeSOCKcluster(5)
			
		
		
		clusterExport(cl, "getCVreplicates.child") 
		#clusterExport(cl, "numsamp")
		#clusterExport(cl, "numreplicates")
		#clusterExport(cl, "min.samp.percent")
		#clusterExport(cl, "impute.bool")
		#clusterExport(cl, "alignment.tool")
		
		cv.res<-parApply(cl,curdata,1,getCVreplicates.child,numsamp=numsamp,numreplicates=numreplicates,min.samp.percent=min.samp.percent,impute.bool=impute.bool,
		missingvalue=missingvalue)
		
	dim(cv.res)=dim(matrix(nrow=7,ncol=numfeats))
	#cv.res<-t(cv.res)
			
		
	stopCluster(cl)
        final_set<-as.data.frame(cv.res)
        rownames(final_set)=NULL
	final_set<-apply(final_set,2,as.numeric)
	final_set<-as.data.frame(t(final_set))

	numsamp=numsamp/numreplicates
	colnames(final_set)<-c("min", "first_quartile", "median", "mean", "third_quartile", "max","numgoodsamples")
	
	#deltappm_res<-apply(curdata_mz_rt_info,1,get_deltappm)
	#delta_cv_range<-as.numeric(final_set$max)-as.numeric(final_set$min)+0.1
	
	#Qscore<-100*((final_set$numgoodsamples)/(delta_cv_range*numsamp*(deltappm_res+0.1)))
	
	
	#Qscore<-100*(final_set$numgoodsamples/(final_set$median*numsamp))

	 mz_min_max<-cbind(curdata_mz_rt_info[,1],curdata_mz_rt_info[,1])	
	if(alignment.tool=="apLCMS"){
		mz_min_max<-cbind(curdata_mz_rt_info[,3],curdata_mz_rt_info[,4])
	}else{
		if(alignment.tool=="XCMS"){
			mz_min_max<-cbind(curdata_mz_rt_info[,2],curdata_mz_rt_info[,3])
		}
	}
	mz_min_max<-as.data.frame(mz_min_max)
	
	deltappm_res<-apply(mz_min_max,1,get_deltappm)
	delta_cv_range<-1 #as.numeric(final_set$max)-as.numeric(final_set$min)+0.1
	
	
	#Qscore<-100*(final_set$numgoodsamples)/((((as.numeric(final_set$median)+1)*numsamp*((10*(deltappm_res+0.1))))+1))
	
	if(is.na(missingvalue)==TRUE){
							num_tot_nonzeros<-apply(curdata,1,function(x){length(which(is.na(x)==FALSE))})
                        
                        
						}else{
							
							num_tot_nonzeros<-apply(curdata,1,function(x){length(which(x>missingvalue))})
                        	
						}
	
	#Qscore<-100*((final_set$numgoodsamples)/(delta_cv_range*as.numeric(final_set$median)*num_tot_zeros*(deltappm_res+0.1)))
	
	#Qscore<-100*((final_set$numgoodsamples)/((as.numeric(final_set$median)*num_tot_zeros*(deltappm_res+0.1))+1)
	
	#Qscore<-100*((final_set$numgoodsamples)/(((as.numeric(final_set$median)+1)*(num_tot_zeros+0.1)*(deltappm_res+0.1))+1))
	
	#Qscore<-100*(final_set$numgoodsamples)/((((as.numeric(final_set$median)+1)*numsamp*(((deltappm_res+0.1))))+1))
	
	#Qscore<-(100*final_set$numgoodsamples)/((((as.numeric(final_set$median)+0.01)*numsamp*(((deltappm_res+0.01))))+1))
	
	#Qscore<-(100*final_set$numgoodsamples)/((as.numeric(final_set$median)+0.01)*num_tot_nonzeros*0.5*(deltappm_res+0.01))
	
	#part of xMSanalyzer v2.0.4
	Qscore<-(100*final_set$numgoodsamples)/(((as.numeric(final_set$median)+0.01+1)*num_tot_nonzeros*(1/numreplicates)*(deltappm_res+0.01+1))+0.01)
	
	termA<-final_set$median+0.01+1
	termB<-num_tot_nonzeros*(1/numreplicates)
	#Qscore<-(100*final_set$numgoodsamples)/(((max(as.numeric(final_set$median),as.numeric(final_set$max),na.rm=TRUE)+0.01+1)*num_tot_nonzeros*(1/numreplicates))+0.01)
	
	#part of xMSanalyzer_v2.0.6
	Qscore<-(100*final_set$numgoodsamples)/((termA*termB)+0.01)
	
	
	
	final_set<-cbind(final_set,Qscore)
	final_set<-as.data.frame(final_set)
	
	#quantcolnames
        return(final_set)
}
