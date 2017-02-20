check.mz.in.replicates.child <-
function(curdata,min.samps,min.reps,num_replicates)
{
        
        numsamp=length(curdata)
        textp1=""
        t=1

        finalmat={}
      
	
                replicate_check=0
		intvec={}
		num.samps.check=0
                for(samp in seq(1,(numsamp),num_replicates))
                {
                        newrow={}
                	intvec={}
		        for(replicate in 1:num_replicates)
			{ 
				i=samp+replicate-1
				intvec=c(intvec,curdata[i])
				
			}
			if(length(which(intvec>0))>=(min.reps))
			{
                                	replicate_check=1
				 	num.samps.check=num.samps.check+1
			}
                
		}
		if(num.samps.check>=min.samps)
		{	
			replicate_check=1
		}
		else
		{
			replicate_check=0
		}
               
	return(replicate_check)

}
