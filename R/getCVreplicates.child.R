getCVreplicates.child <-
function(curdata,numsamp,numreplicates,min.samp.percent,impute.bool,missingvalue)
{
	
			newrow={}
			#numsamp=length(curdata)
			goodsamps<-0
			for(samp in seq(1,numsamp,numreplicates))
			{
				i=samp
				j=i+numreplicates-1
				
				curdata_int=curdata[c(i:j)]
				
				if(is.na(missingvalue)==TRUE){
							
                        check_zeros=which(is.na(curdata_int)==TRUE)
                        
						}else{
							check_zeros=which(curdata_int==missingvalue)
								
						}
				#check_zeros=which(curdata_int==0)
				
				na_thresh=round(min.samp.percent*numreplicates)
				
			
					if(length(check_zeros)>=na_thresh)
					{
						
							cvval<-NA
						
						#newrow<-cbind(newrow,cvval)
					
					}
					else
					{
						#temporarily replace the missing intensities, set to 0 in apLCMS,
						#with mean intensity value of the corresponding replicates (with non-zero values)
						if(length(check_zeros)>0){
						if(impute.bool==TRUE){
							curdata_int[check_zeros]=mean(t(curdata_int[-c(check_zeros)]))
						}
						}
						sdval<-sd(curdata_int)
						meanval<-mean(t(curdata_int))
						cvval<-100*(sdval/meanval)
						newrow<-cbind(newrow,cvval)
						
						goodsamps<-goodsamps+1
					}
				
			}
			if(length(newrow)>0)
			{
				na_ind=which(is.na(newrow)==TRUE)
			
				if(length(na_ind)>0)
				{
					sumrow=summary(as.vector(newrow[-na_ind]))
				}
				else
				{
					sumrow=summary(as.vector(newrow))
				}
			}
			else{sumrow<-{}}
			if(length(sumrow)<6)
			{
				for(i in 1:6)
				{
				 sumrow[i]=200
				}
			}
			finalmat<-{}
			
			finalmat<-rbind(finalmat, c(unlist(sumrow),goodsamps))
			return(finalmat)

}
