merge.Results.child.ttest <-
function(dataA, max.mz.diff=15, max.rt.diff=300, merge.eval.pvalue=0.05,alignment.tool="apLCMS", mult.test.cor=FALSE,mergecorthresh=0.7,col.rm.index=NA,missingvalue=0)
{	
        diff_mz_num=1
        unique_mz={}
        ppm_v={}
        rt_v={}
        cnames=colnames(dataA)
        dataA<-as.data.frame(dataA)
        if(nrow(dataA)<2){
        	
        	return(dataA)
        }
        
        if(is.na(col.rm.index)==FALSE){
       	 dataA<-as.data.frame(dataA[,-c(col.rm.index)]) 
		}else{
			 dataA<-as.data.frame(dataA)
			}
			
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

		if(is.na(missingvalue)==FALSE){
		
			dataA<-replace(as.matrix(dataA),which(dataA==missingvalue),NA)
			
		}
		
		dataA<-as.data.frame(dataA)
  
  compare_intensities_ttest<-function(other_feats,y,merge.eval.pvalue, mergecorthresh=0.6){
		
											ttest_sum=apply(other_feats,1,function(x)
				                                                        {
				                                                        	
												diff_rt<-abs(x[1]-y[1])
												
												x<-x[-c(1)]
												y<-y[-c(1)]
												x<-as.matrix(x)
				                                y<-as.matrix(y)
				                                
				                                
												yind<-which(is.na(y)==TRUE)
												xind<-which(is.na(x)==TRUE)
												
												naind<-c(yind,xind)
												
				                                if(dim(x)[1]>dim(y)[1])
												{
													x<-t(as.numeric(x))
												}
												
												if(dim(y)[1]>dim(x)[1])
												{
													y<-t(as.numeric(y))
												}
												
												if(length(naind)>0){
												 x<-x[-naind]
												 y<-y[-naind]
				                                                                 
												 }
												
												 if(length(x)>1 & length(y)>1)
												 {
													num_com<-length(which(abs(x-y)<0.5))
													
													if(num_com<1)
													{
														
													#	print(x)
													#	print(y)
													
													mean_x<-mean((x+1),na.rm=TRUE)
													mean_y<-mean((y+1),na.rm=TRUE)
													
												#	print(mean_x)
												#	print(mean_y)
													log2fc<-abs(log2(mean_x/mean_y))
													
												#	print(log2fc)
													
													if(log2fc>0.1)
													{
													
													
													ttest_res=try(t.test(x,y,paired=TRUE,na.omit=TRUE),silent=TRUE)
													
													#print(ttest_res)
													if (is(ttest_res, "try-error"))
													{
														ttest_pval<-0
														ttest_pval2<-0

													}else
													{
														
														ttest_pval=ttest_res$p.value
														
														if(is.na(ttest_pval)==TRUE){
															ttest_pval=1
														}
														
														ttest_res2=try(t.test(x,y,paired=FALSE,na.omit=TRUE),silent=TRUE)
														ttest_pval2=ttest_res2$p.value
														if(is.na(ttest_pval2)==TRUE){
															ttest_pval2=1
														}
														
													}
													
													
													
													
													#if(ttest_pval<merge.eval.pvalue)
													
													
													#r6<-kruskal.test(list(x,y))
													
													#r10<-var.test(x,y)
													
													ttest_pval<-max(ttest_pval,ttest_pval2,na.rm=TRUE)
													
													if(length(x)>2 & length(y)>2)
													{
																cortest_res=try(cor.test(as.numeric(x),as.numeric(y),method="spearman",use="pairwise.complete.obs"),silent=TRUE)
																if (is(cortest_res, "try-error")){
																	cortest_pval<-1
																	cortest_est<-0
																}else{
																	cortest_pval=cortest_res$p.value
																	cortest_est=cortest_res$estimate
																	
																	if(cortest_est>=mergecorthresh){
																		
																		ttest_pval<-1
																		cortest_est<-1
																		cortest_pval<-0
																		
																	}
																	
																	if(cortest_est<mergecorthresh){
																		cortest_pval=1
																		
																					cortest_res=try(cor.test(as.numeric(x),as.numeric(y),method="pearson",use="pairwise.complete.obs"),silent=TRUE)
																				if (is(cortest_res, "try-error")){
																					
																				cortest_pval=1
																				cortest_est<-0
																				}else{
																					
																					cortest_pval=cortest_res$p.value
																					cortest_est=cortest_res$estimate
																					
																					if(cortest_est>=mergecorthresh){
																		
																		ttest_pval<-1
																		cortest_est<-1
																		cortest_pval<-0
																		
																	}
																					
																					if(cortest_est<mergecorthresh){
																						cortest_pval=1
																						cortest_est<-0
																					}
																					
																					}
																	}
																}
													
													
													
													if(diff_rt>30)
													{

													if(cortest_est>=mergecorthresh & ttest_pval>=merge.eval.pvalue)
                                                    {
													
													
															ttest_pval=1
													
													}else{
														ttest_pval=0
														}
													}else{
														
														if(cortest_est>=mergecorthresh | ttest_pval>=merge.eval.pvalue)
														{
													
																#print("here")
																ttest_pval=1
														}
														
														
														
														if(diff_rt<2)
														{
															ttest_pval=1
															
														}
															
														
														}
													
													}
													else{
														
														ttest_pval<-ttest_pval
													}
													#ttest_pval<-max(c(ttest_pval,r6$p.value,r10$p.value),na.rm=TRUE)
													
													
													}
													else
													{
													
														ttest_pval=1
													}
													
													
													}
													else
													{
													
													ttest_pval=1
													}											
													
													
				                                  }
				                                 else
				                                 {
				                                 	
				                                 	#if(((x/y)>2) | (x/y)<0.5){
				                                 		
				                                 	#	ttest_pval<-0
				                                 	#}else{
														ttest_pval=1
													#}
												}
												
				                                                                return(ttest_pval)
				                                                        })
return(ttest_sum)
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
          
    fname<-paste("mz_groups",max.mz.diff,".Rda",sep="")
     
    #save(mz_groups,file=fname)

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

	
	
        #Step 2 Sub-group features from Step 1 by Retention time
        #find the features with RT values within the defined range as compared to the query feature
         diffmat={} #length(mz_groups)
       # for(j in 1:length(mz_groups))
       diffmatres<-lapply(1:length(mz_groups),function(j)
        {
        	fname<-paste("mz_groups",max.mz.diff,"_",j,".Rda",sep="")
     
    			
               temp_diff={}
                tempdata=mz_groups[[j]][[1]]
              # save(tempdata,file=fname)
		
                if(dim(tempdata)[1]>1)
                {
				
                #print(tempdata[,1:2])
				
				
                                tempdata<-tempdata[order(tempdata$Qscore,decreasing=TRUE),]
				
				#rem_ind<-apply(tempdata,1,function(x){as.numeric(x[(dim(tempdata)[2]-7)])==as.numeric(x[(dim(tempdata)[2]-2)])})
				
				rem_ind<-{}
				
				rem_ind<-which(rem_ind==TRUE)
								if(length(rem_ind)>0){
								
									
									 tempdata<-tempdata[-rem_ind,]
								}
				allrownames<-rownames(tempdata)
					
				dup_list<-{}
				temp_diff={}
				dup_ind<-{}
				curdupnames<-{}
				diffmat<-{}
                                for(d in 1:dim(tempdata)[1])
                                {
                                		
										if(d%in%dup_ind==FALSE & d%in%curdupnames==FALSE)
										{
										rowname<-rownames(tempdata[d,])                                      
					                     cur_feat=tempdata[d,]
										
										 other_feats=tempdata
									
										getbind_rtsame<-which(abs(other_feats$time-cur_feat$time)<=max.rt.diff)
					                                        commat={}
										 same_feat_ind={}
									
										 
					                    if(length(getbind_rtsame)>1)
                                        {
                                                other_feats<-other_feats[getbind_rtsame,]
		
												y=cur_feat[c(2,sample.col.start:(dim(tempdata)[2]-8))]
						
												diff_rt<-abs(other_feats[,2]-y[2])

												ttest_sum<-compare_intensities_ttest(other_feats[,c(2,sample.col.start:(dim(tempdata)[2]-8))],y,merge.eval.pvalue,mergecorthresh)
											
					                            same_feat_ind=which(ttest_sum>(merge.eval.pvalue/num_comparisons))
					                            same_feat_ind=c(same_feat_ind,which(is.na(ttest_sum)))
												dup_list<-c(dup_list,names(same_feat_ind))
												curnames<-names(same_feat_ind)
												curdupnames<-which(allrownames%in%curnames)
												
												curnames<-curnames[-which(curnames==rowname)]
												
												
													
												if(length(which(curdupnames%in%dup_ind==TRUE))<1)
													{
																dup_ind<-c(dup_ind,curdupnames)
														
					                                                        if(length(same_feat_ind)>1)
					                                                        {
					                                                                                        same_ind<-as.numeric(same_feat_ind)						
													other_feats<-as.data.frame(other_feats)
									
					                                                                commat=other_feats[c(same_ind),]
					                                                               maxint<-apply(commat,1,function(x){max(x,na.rm=TRUE)})
					                                                                
					                                                                maxint<-log((maxint+1),10)
					                                                                maxint<-maxint+0.001
					                                                                maxint<-maxint/max(maxint,na.rm=TRUE)
					                                                                
						                                           commat_qscore<-as.numeric(commat$Qscore)*maxint #cv_range/as.numeric(commat$numgoodsamples)*as.numeric(deltappm_res)
																	#inverse of the actual qscore
													best_level_index=which(as.numeric(commat_qscore)==max(as.numeric(commat_qscore),na.rm=TRUE))
					                                                                if(length(best_level_index)>0)
					                                                                {
					                                                                        best_level_index=best_level_index[1]
					                                                                }
					                                                                best_level_data=commat[best_level_index,]
					                                                        }else
					                                                        {
																					best_level_index=d
					                                                                best_level_data=cur_feat
					
					                                                        }
												
														
														
					                                                                       diffmat=rbind(diffmat, best_level_data)
					                                                                        #diffmat<-as.data.frame(best_level_data)

					                                                                        temp_diff=rbind(temp_diff,best_level_data)
														
																							dup_ind<-c(dup_ind,best_level_index)
														
													} else{
                                											if(d%in%dup_ind==FALSE){
				                                               					 diffmat<-rbind(diffmat,cur_feat)
				                                               					 #diffmat<-as.data.frame(cur_feat)

				                                              					  temp_diff=rbind(temp_diff,tempdata[d,])
																				  dup_ind<-c(dup_ind,d)
				                                              					  }
                                				}
												
					
                                }
                                else{
                                			if(d%in%dup_ind==FALSE){
				                                               					 diffmat<-rbind(diffmat,cur_feat)
				                                               					 #diffmat<-as.data.frame(cur_feat)

				                                              					  temp_diff=rbind(temp_diff,tempdata[d,])
																				  dup_ind<-c(dup_ind,d)
				                                              					  }
                                	}
                               }
                               
							}#for loop
                	
                }
                else{
                	cur_feat<-as.matrix(mz_groups[[j]][[1]])
                	
                        
                                               					 diffmat<-rbind(diffmat,cur_feat)
                                               					 # diffmat<-as.data.frame(cur_feat)
                                              					 temp_diff=rbind(temp_diff,tempdata)
                        	}
                        
                        
                        
                
                diffmat<-unique(diffmat)
                 best_level_data={}

		return(diffmat)
        }
        )
        
        diffmat<-do.call(rbind,diffmatres)
        diffmat=unique(diffmat)
        return(diffmat)
        
}
