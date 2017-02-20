getPID <-
function(curdata, alignment.tool,missingvalue, numreplicates)
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
        numfeats=dim(curdata)[1]
        numsamp=dim(curdata)[2]
        
        maxint<-apply(curdata,1,function(x){max(x,na.rm=TRUE)})
        maxint<-log(maxint,10)
        
       # maxint<-(maxint)/(max(maxint,na.rm=TRUE))
                
        rnames<-colnames(curdata)
        rnames<-gsub(".cdf", "", rnames, ignore.case=TRUE)
        quantcolnames=c("min", "first_quartile", "median", "mean", "third_quartile", "max")
	
        resvec_1<-lapply(1:numfeats,function(r)
        {
                newrow={}
                finalmat={}
		no_value=0
		goodsamps<-0
                for(samp in seq(1,(numsamp),2))
                {
                        i=samp
                        j=i+1
                        int1=curdata[r,i]
                        int2=curdata[r,j]
			
			curdata_int=curdata[r,c(i:j)]
						if(is.na(missingvalue)==TRUE){
							
                        check_zeros=which(is.na(curdata_int)==TRUE)
                        
						}else{
							check_zeros=which(curdata_int==missingvalue)
								
						}
                        
			 			if(length(check_zeros)>0)
                        {
                                                 
                                replicate_diff<-NA
								no_value<-no_value+1
                        }
                        else
                        {
                                #calculate PID
                                replicate_diff<-100*(abs(int1-int2)/mean(c(int1,int2)))
								goodsamps<-goodsamps+1
                        }
                        newrow<-cbind(newrow,replicate_diff)
                }

                #get indices of the PIDs that are NA
                na_ind=which(is.na(newrow)==TRUE)

                #get quantile summary of the percent intensity difference (PID) vector
                #using only the non-NA values
                if(length(na_ind)>0)
                {
                        sumrow=summary(as.vector(newrow[-na_ind]))
                }
                else
                {
                        sumrow=summary(as.vector(newrow))
                }

                #if quantile values are set  to NA
                if(length(sumrow)<6)
                {
                        for(i in 1:6)
                        {
                                sumrow[i]=200
                        }
                }
                names(sumrow)=quantcolnames
		#tempres<-unlist(sumrow)
		#tempres<-cbind(tempres,goodsamps)
                #finalmat<-rbind(finalmat, tempres)
		
		finalmat<-rbind(finalmat, c(unlist(sumrow),goodsamps))
                return(finalmat)
        })
        final_set={}
        for(i in 1:length(resvec_1)){
                if(length(resvec_1[[i]])>1)
                {
                        final_set<-rbind(final_set,resvec_1[[i]])
                }
        }

        final_set<-as.data.frame(final_set)
        rownames(final_set)=NULL
	    final_set<-apply(final_set,2,as.numeric)
		final_set<-as.data.frame(final_set)
		
		colnames(final_set)<-c("min", "first_quartile", "median", "mean", "third_quartile", "max", "numgoodsamples")
		numsamp<-numsamp/2
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
	delta_cv_range<-as.numeric(final_set$max)-as.numeric(final_set$min)+0.1
	
	#Qscore<-100*((final_set$numgoodsamples)/(delta_cv_range*as.numeric(final_set$median)*numsamp*(deltappm_res+0.1)))
	
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
	
	#Qscore<-(100*final_set$numgoodsamples)/(((as.numeric(final_set$median)+0.01)*num_tot_nonzeros*0.5*(deltappm_res+0.01))+0.01)
	
	#Qscore<-Qscore*(final_set$numgoodsamples/numsamp)
	
	
	#part of xMSanalyzer_v2.0.4
	#Qscore<-(100*final_set$numgoodsamples)/(((as.numeric(final_set$median)+0.01+1)*num_tot_nonzeros*0.5*(1+deltappm_res+0.01))+0.01)
	
	#part of xMSanalyzer_v2.0.5
	#Qscore<-(100*final_set$numgoodsamples)/(((max(as.numeric(final_set$median),as.numeric(final_set$max),na.rm=TRUE)+0.01+1)*num_tot_nonzeros*0.5)+0.01)
	
	#part of xMSanalyzer_v2.0.6
	termA<-final_set$median+0.01+1
        termB<-num_tot_nonzeros*(1/numreplicates)

	Qscore<-(100*final_set$numgoodsamples)/((termA*termB)+0.01)

	final_set<-cbind(final_set,Qscore)
	final_set<-as.data.frame(final_set)
        return(final_set)
}
