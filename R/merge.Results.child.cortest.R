merge.Results.child.cortest <-
function(dataA, max.mz.diff=15, max.rt.diff=300, merge.eval.pvalue=0.05,alignment.tool="apLCMS", mult.test.cor=FALSE)
{
    
        diff_mz_num=1
        unique_mz={}
        ppm_v={}
        rt_v={}
        cnames=colnames(dataA)
        dataA<-as.data.frame(dataA)
	
		#dataA<-as.data.frame(dataA,2,as.numeric)
         if(alignment.tool=="apLCMS")
        {
              sample.col.start=6
        }
        else
        {
              if(alignment.tool=="XCMS")
              {
                    sample.col.start=9
                    cnames[1]="mz"

                    cnames[4]="time"
                    colnames(dataA)=cnames
                    
              }
              else
              {
                    stop(paste("Invalid value for alignment.tool. Please use either \"apLCMS\" or \"XCMS\"", sep=""))
              }
        }

compare_intensities_cortest<-function(other_feats,y){
											cortest_sum=apply(other_feats,1,function(x)
				                                                        {
				                                                        	x<-as.matrix(x)
				                                                        	y<-as.matrix(y)
												yind<-which(y==0)
												xind<-which(x==0)
												naind<-c(yind,xind)
												
												
				                                                                if(dim(x)[1]>dim(y)[1])
												{
													x<-t(x)
												}
												if(dim(y)[1]>dim(x)[1])
												{
													y<-t(y)
												}
				                                                                      if(length(naind)>0){
														 x<-x[-naind]
														 y<-y[-naind]
				                                                                 
													}
												
													#if(max(x)!=0 & max(y)!=0)
													if(length(x)>2 & length(y)>2)
													{
														
														
														cortest_res=try(cor.test(x,y),silent=TRUE)
																if (is(cortest_res, "try-error")){
																	cortest_pval<-1
																	cortest_est<-0
																}else{
																cortest_pval=cortest_res$p.value
																cortest_est=cortest_res$estimate
																if(cortest_est<0){
																	cortest_pval=1
																}
																}
													}else{
													cortest_res=try(cor.test(x,y),silent=TRUE)
																if (is(cortest_res, "try-error")){
																	cortest_pval<-1
																	cortest_est<-0
																}else{
																cortest_pval=cortest_res$p.value
																cortest_est=cortest_res$estimate
																if(cortest_est<0){
																	cortest_pval=1
																}
																}
													}
												
												
				                                                                return(cortest_pval)
				                                                        })
return(cortest_sum)
}


        #Step 1 Group features by m/z
        mz_groups<-lapply(1:dim(dataA)[1],function(j){
                                commat={}
                                diffmz=new("list")
                                ppmb=(max.mz.diff)*(dataA$mz[j]/1000000)
                                getbind_same<-which(abs(dataA$mz-dataA$mz[j])<=ppmb)
                                diffmz[[diff_mz_num]]=getbind_same #dataA[getbind_same,]
                                diff_mz_num=diff_mz_num+1
                                return(diffmz)
          })
	  
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
	
	diff_mz_num=1
	mz_groups<-lapply(1:length(mz_groups),function(j){
                                commat={}
                                diffmz=new("list")
                                getbind_same=mz_groups[[j]][[1]]
                                diffmz[[diff_mz_num]]=dataA[getbind_same,]
                                diff_mz_num=diff_mz_num+1
                                return(diffmz)
          })
	  
	if(mult.test.cor==TRUE){
		mz_group_size<-lapply(1:length(mz_groups),function(j){
		num_rows<-dim(mz_groups[[j]][[1]])
		n=num_rows[[1]][[1]]
		num_comp=(n*(n-1))/2
		})

	
		num_comparisons<-sum(unlist(mz_group_size))
	}else{
		num_comparisons=1
	}
	print("num comparisons")
	print(num_comparisons)
	
        #Step 2 Sub-group features from Step 1 by Retention time
        #find the features with RT values within the defined range as compared to the query feature
        diffmat={}
        for(j in 1:length(mz_groups))
        {
                temp_diff={}
                tempdata=mz_groups[[j]][[1]]
		
                if(dim(tempdata)[1]>1)
                {
		
                                tempdata<-tempdata[order(tempdata$median),]
				
				rem_ind<-apply(tempdata,1,function(x){as.numeric(x[(dim(tempdata)[2]-6)])==as.numeric(x[(dim(tempdata)[2]-1)])})
				
				rem_ind<-which(rem_ind==TRUE)
								if(length(rem_ind)>0){
								
									
									 tempdata<-tempdata[-rem_ind,]
								}
				if(dim(tempdata)[1]>0){
				dup_list<-{}
				temp_diff={}
                                for(d in 1:dim(tempdata)[1])
                                {
					rowname<-rownames(tempdata[d,])
					print(rowname)
					print(tempdata[,1:5])
					print(length(tempdata))
					if((rowname%in%dup_list)==FALSE)
					{
                                        #tempdata=mz_groups[[j]][[1]]
                                        cur_feat=tempdata[d,]
                                        other_feats=tempdata
                                        getbind_rtsame<-which(abs(other_feats$time-cur_feat$time)<=max.rt.diff)
                                        commat={}
					 same_feat_ind={}
                                        if(length(getbind_rtsame)>1)
                                        {
                                                other_feats<-other_feats[getbind_rtsame,]
		
					y=cur_feat[sample.col.start:(dim(tempdata)[2]-7)]
		
										 
										ttest_sum<-compare_intensities_cortest(other_feats[,sample.col.start:(dim(tempdata)[2]-7)],y)

                                                        same_feat_ind=which(ttest_sum<(merge.eval.pvalue/num_comparisons))
                                                        same_feat_ind=c(same_feat_ind,which(is.na(ttest_sum)))
							print(ttest_sum)
							dup_list<-c(dup_list,names(same_feat_ind))
									 				
                                                        if(length(same_feat_ind)>0)
                                                        {
                                                                commat=other_feats[same_feat_ind,]
								#commat_qscore<-as.numeric(commat$median)/apply(commat[,sample.col.start:(dim(tempdata)[2]-7)],1,countpeaks)
                                                                #best_level_index=which(as.numeric(commat$median)==min(as.numeric(commat$median)))
								
								 commat_qscore<-as.numeric(commat$median)/as.numeric(commat$numgoodsamples)
								 
								best_level_index=which(as.numeric(commat_qscore)==min(as.numeric(commat_qscore)))
                                                                if(length(best_level_index)>0)
                                                                {
                                                                        best_level_index=best_level_index[1]
                                                                }
                                                                best_level_data=commat[best_level_index,]
                                                        }else
                                                        {
                                                                best_level_data=cur_feat

                                                        }
                          
                                                                        diffmat=rbind(diffmat, best_level_data)
                                                                        temp_diff=rbind(temp_diff,best_level_data)
					
                                        }
                                        else
                                        {
          

                                               					 diffmat<-rbind(diffmat,cur_feat)
                                              					  temp_diff=rbind(temp_diff,tempdata[d,])
                                              					  
                                        		
                                        		
                                        }
					}
                                }
				}
                }
                else
                {
                	cur_feat<-as.matrix(mz_groups[[j]][[1]])
                	
                        									#tempdata=mz_groups[[j]][[1]]
                         if(tempdata$min!=tempdata$max)
			 {
                                               					 diffmat<-rbind(diffmat,cur_feat)
                                              					 temp_diff=rbind(temp_diff,tempdata)
                        	}
                        
                        
                        
                }
                diffmat<-unique(diffmat)
                 best_level_data={}
                

        }
        diffmat=unique(diffmat)
        return(diffmat)
}
