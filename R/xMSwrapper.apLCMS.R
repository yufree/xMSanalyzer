xMSwrapper.apLCMS <-
function(cdfloc, apLCMS.outloc, xMSanalyzer.outloc, min.run.list=c(4,3), min.pres.list=c(0.5,0.8), minexp.pct=0.1, mztol=NA, alignmztol=NA, 
alignchrtol=NA,numnodes=NA,run.order.file=NA, apLCMSmode="untargeted", known_table=NA, match_tol_ppm=5, max.mz.diff=15,max.rt.diff=300,merge.eval.pvalue=0.2,mergecorthresh=0.7,deltamzminmax.tol=10, subs=NA ,num_replicates=3, 
mz.tolerance.dbmatch=10, adduct.list=c("M+H"), samp.filt.thresh=0.70,feat.filt.thresh=50,cormethod="pearson", mult.test.cor=TRUE,
missingvalue=0,ignore.missing=TRUE,filepattern=".cdf",sample_info_file=NA,refMZ=NA,refMZ.mz.diff=10,refMZ.time.diff=NA,void.vol.timethresh=30,replacezeroswithNA=TRUE,
scoreweight=30,charge_type="pos",syssleep=0.5)
{

	suppressWarnings(suppressWarnings(sink(file=NULL)))
	x<-date()
	x<-strsplit(x,split=" ")
	
	targeted_feat_raw<-{}
	
	x1<-unlist(x)
	#fname<-paste(x1[2:5],collapse="_")
	
	#fname<-paste(x1[2:3],x1[5],x1[4],collapse="_")
	
	fname<-paste(x1[2:3],collapse="")
	fname<-paste(fname,x1[5],sep="")
	x1[4]<-gsub(x1[4],pattern=":",replacement="_")
	fname<-paste(fname,x1[4],sep="_")
	
	
	fname<-paste(xMSanalyzer.outloc,"/Log_",fname,".txt",sep="")
	print(paste("Program is running. Please check the logfile for runtime status: ",fname,sep=""))

	 data_rpd_all=new("list")
	data_rpd=new("list")
	union_list=new("list")
	feateval_list<-new("list")
	rsqres_list<-new("list")
	annot.res<-{}
	annotres_list<-new("list")
	minexp<-2
	 if(is.na(sample_info_file)==FALSE)
    {
    	
    	
		sampleid_mapping<-read.table(sample_info_file,sep="\t",header=TRUE)

		if(is.na(cdfloc)==FALSE){
			l1<-list.files(cdfloc,filepattern)

			if(is.na(minexp.pct)==FALSE){	
			minexp<-round(minexp.pct*length(l1))
			}else{
			minexp<-2

			}
		}else{
				l1<-rep(1,num_replicates)	
			}
		if(length(l1)!=dim(sampleid_mapping)[1] & (is.na(cdfloc)==FALSE))
		{
			num_mis_files<-dim(sampleid_mapping)[1]-length(l1)
			stop(paste("ERROR: Only ",length(l1)," spectral files were found. ",num_mis_files," files are missing.",sep=""))
		
		}
	}else{
	
	
		if(is.na(cdfloc)==FALSE){
			l1<-list.files(cdfloc,filepattern)
			
			 if(is.na(minexp.pct)==FALSE){
                        minexp<-round(minexp.pct*length(l1))
                        }else{
                        minexp<-2

                        }
			
		
		}else{
			l1<-rep(1,num_replicates)	
		}
	}
if(length(l1)%%num_replicates>0)
{stop(paste("ERROR: Not all samples have ",num_replicates," replicates.",sep=""))
}


dir.create(apLCMS.outloc,showWarnings=FALSE)
	dir.create(xMSanalyzer.outloc,showWarnings=FALSE)
	
	
	sink(fname)
	print(sessionInfo())
	    if(is.na(refMZ)==FALSE){
			stddata<-read.table(refMZ,sep="\t",header=TRUE)
			print(refMZ)
			print(head(stddata))
			
		}else{
			if(charge_type=="pos"){
			data(example_target_list_pos)
			stddata<-example_target_list_pos
			}else{
				
				if(charge_type=="neg"){
						data(example_target_list_neg)
						stddata<-example_target_list_neg
					}else{
						stop("Invalid option. \'charge_type\' should be \'pos\' or \'neg\'.")
						}
				}
		}
        ############################################
        #1) Align profiles using the cdf.to.ftr wrapper function in apLCMS
        if(is.na(apLCMS.outloc)==TRUE)
	{
		stop("Undefined value for parameter, apLCMS.outloc. Please define the apLCMS output location.")

	}
	 if(is.na(xMSanalyzer.outloc)==TRUE)
        {
                stop("Undefined value for parameter, xMSanalyzer.outloc. Please define the xMSanalyzer output location.")

        }

	if(is.na(cdfloc)==FALSE)
        {
                setwd(cdfloc)
                if(is.na(apLCMS.outloc)==FALSE)
                {
                      
			data_rpd_all<-apLCMS.align(cdfloc, apLCMS.outloc,min.run.list, min.pres.list, minexp=minexp, mztol=mztol, alignmztol=alignmztol, alignchrtol=alignchrtol, 
			numnodes,run.order.file,subs,filepattern=filepattern,apLCMSmode=apLCMSmode, known_table=known_table, match_tol_ppm=match_tol_ppm,
			refMZ.mz.diff,refMZ.time.diff,target.mz.list = stddata)
    
                }
                else
                {
                        stop("Undefined value for parameter, apLCMS.outloc. Please define the output location.")
                }
        }
	
		setwd(apLCMS.outloc)
                alignmentresults<-list.files(apLCMS.outloc, "*.txt")
		print("Files found in apLCMS output location:")
		print(alignmentresults)
		for(i in 1:length(alignmentresults)){
		
			data_rpd_all[[i]]<-read.table(paste(apLCMS.outloc,"/",alignmentresults[i],sep=""),header=TRUE)
			#data_rpd_all[[i]]<-data_rpd_all[[i]]
			
			#if(is.na(cdfloc)==TRUE)
			{
			l1<-dim(data_rpd_all[[i]])[2]-4
	
			}
			
			
			numfeats<-apply(data_rpd_all[[i]][,-c(1:4)],1,function(x){
					return(countpeaks(x, missingvalue))
					})
					numfeats<-unlist(numfeats)
					 if(is.na(minexp.pct)==FALSE)
					{
					minexp<-minexp.pct*l1
					
					if(length(which(numfeats>=minexp))>0){
					print(paste("Removing ",length(which(numfeats<minexp))," features based on missing value criteria",sep=""))
					data_rpd_all[[i]]<-data_rpd_all[[i]][which(numfeats>=minexp),]
					}else{
						stop(paste("No features have non-missing value in ",minexp, " samples",sep=""))
						}
					
					
					}
			
		}
	
	
	
	#subdir1<-paste(xMSanalyzer.outloc,"/Quality_assessment_files",sep="")
	#subdir2<-paste(xMSanalyzer.outloc,"/apLCMS_filtered_data",sep="")
	#subdir3<-paste(xMSanalyzer.outloc,"/apLCMS_with_xMSanalyzer_merged_data",sep="")
	
	#dir.create(subdir1,showWarnings=FALSE)
	#dir.create(subdir2,showWarnings=FALSE)
	#dir.create(subdir3,showWarnings=FALSE)
        
    subdir1<-paste(xMSanalyzer.outloc,"/Stage1",sep="")  #QC individual parameter settings
	subdir2<-paste(xMSanalyzer.outloc,"/Stage2",sep="")  #Data filtering
	subdir3<-paste(xMSanalyzer.outloc,"/Stage3a",sep="")  #Data merger/parameter optimization
	subdir3b<-paste(xMSanalyzer.outloc,"/Stage3b",sep="") 
	subdir4a<-paste(xMSanalyzer.outloc,"/Stage4a",sep="")	 #Raw QC: batch effect eval, TIC, etc
	
	
	dir.create(subdir1,showWarnings=FALSE)
	dir.create(subdir2,showWarnings=FALSE)
	dir.create(subdir3,showWarnings=FALSE)
	dir.create(subdir3b,showWarnings=FALSE)
	dir.create(subdir4a,showWarnings=FALSE)
	
	  if(is.na(sample_info_file)==FALSE)
    {
		subdir4b<-paste(xMSanalyzer.outloc,"/Stage4b",sep="")	 #Batch-effect corrected QC: batch effect eval, TIC, etc
		dir.create(subdir4b,showWarnings=FALSE)
		
	}
	
	if(is.na(adduct.list)==FALSE){
	subdir5<-paste(xMSanalyzer.outloc,"/Stage5",sep="")	 #Putative unprocessed annotations;
	

	dir.create(subdir5,showWarnings=FALSE)
	}
	
	bestscore<-(-1000000)
	
                #stop("Undefined value for parameter, cdfloc. Please enter path of the folder where the CDF files to be processed are located.")
                #change location to the output folder
                setwd(apLCMS.outloc)
                alignmentresults<-list.files(apLCMS.outloc, "*.txt")
                if(length(alignmentresults)>0)
                {
                          curdata_dim={}
                         
			  
			  parent_bad_list<-{}
			  
                          if(num_replicates==2)
                          {
                                  fileroot="_PID"
                          }
                          else
                          {
                                  if(num_replicates>2)
                                  {
                                          fileroot="_CV"
                                  }
                                  else
                                  {
                                          fileroot=""
                        		stop("Need at least 2 technical replicates per sample.")  
			        }
                          }
			  
			  pfilenames<-{}
			  
			   print("*******xMSanalyzer Stage 1: QC evaluation of invidual parameters*******")
			   
			   cat("\n")
                          for(i in 1:length(data_rpd_all))
                          {
                          	
				  print(paste("*Evaluating apLCMS results from parameter setting ",i,"*",sep=""))				  
                                  ############################################
                                  #2)Calculate pairwise correlation coefficients
                                  file_name=sapply(strsplit(alignmentresults[i],".txt"),head)
				  
				  print(paste("*File name: ",alignmentresults[i],"*",sep=""))	
				  
				  crow<-c(paste("p",i,sep=""),file_name)
				  pfilenames<-rbind(pfilenames,crow)

                                  curdata=data_rpd_all[[i]]
                                  #############################################
                                  ############################################
                                  #3) Calculate Percent Intensity Difference

                               if(num_replicates>1)
                                  {
				  
                  
                                feat.eval.result=evaluate.Features(curdata, numreplicates=num_replicates,min.samp.percent=0.6,alignment.tool="apLCMS",impute.bool=TRUE)
                                
                                cnames=colnames(feat.eval.result)
                                feat.eval.result<-apply(feat.eval.result,2,as.numeric)
                                feat.eval.result<-as.data.frame(feat.eval.result)
                                feat.eval.result.mat=cbind(curdata[,c(1:4)],feat.eval.result)
                                feat.eval.outfile=paste(subdir1,"/",file_name,fileroot,"featureassessment.txt",sep="")
                                #write results
                                write.table(feat.eval.result.mat, feat.eval.outfile,sep="\t", row.names=FALSE)
                                
                                curdata<-curdata[which(as.numeric(feat.eval.result$median)<=feat.filt.thresh),]
                                
                                feat.eval.result.mat<-feat.eval.result.mat[which(as.numeric(feat.eval.result$median)<=feat.filt.thresh),]
                                
                                #curdata<-cbind(curdata,feat.eval.result[which(as.numeric(feat.eval.result$median)<=feat.filt.thresh),])
                                
                                curdata<-as.data.frame(curdata)
                                curdata<-replace(as.matrix(curdata),which(is.na(curdata)==TRUE),0)
                                
                                if(is.na(deltamzminmax.tol)==FALSE){
                                    print("filtering by delta m/z")
                                    mz_min_max<-cbind(curdata[,3],curdata[,4])
                                    mz_min_max<-as.data.frame(mz_min_max)
                                    
                                    deltappm_res<-apply(mz_min_max,1,get_deltappm)
                                    
                                    curdata<-curdata[which(as.numeric(deltappm_res)<=deltamzminmax.tol),]
                                    feat.eval.result.mat<-feat.eval.result.mat[which(as.numeric(deltappm_res)<=deltamzminmax.tol),]
                                }
                                
                                feateval_list[[i]]<-feat.eval.result.mat
                                data_rpd[[i]]<-curdata
								
								  if(num_replicates>1)
								  {
									  print(paste("**calculating pairwise ",cormethod," correlation**",sep=""))

									  
									  
											rsqres_list<-evaluate.Samples(curdata, num_replicates, alignment.tool="apLCMS", cormethod,missingvalue,ignore.missing)

											rsqres<-as.data.frame(rsqres_list$cor.matrix)
											
											curdata<-as.data.frame(rsqres_list$feature.table)
											snames<-colnames(curdata[,-c(1:4)])
											snames_1<-snames[seq(1,length(snames),num_replicates)]
											rownames(rsqres)<-snames_1
											pcor_outfile=paste(subdir1,"/",file_name,"_sampleassessment_usinggoodfeatures.txt",sep="")
											write.table(rsqres, pcor_outfile,sep="\t",row.names=TRUE)
											rsqres_list[[i]]<-rsqres
								  }
								  else
								  {
									  print("**skipping sample evaluataion as only one replicate is present**")
								  }
									
								if(num_replicates>2)
								{
									bad_samples<-which(rsqres$meanCorrelation<samp.filt.thresh)
								}else
								{
									bad_samples<-which(rsqres$Correlation<samp.filt.thresh)
								}
								
								if(length(bad_samples)>0){
									bad_sample_names<-snames_1[bad_samples]
									
									feat.eval.outfile=paste(subdir1,"/",file_name,"_badsamples_at_cor",samp.filt.thresh,".txt",sep="")
									bad_sample_names<-as.data.frame(bad_sample_names)
									colnames(bad_sample_names)<-paste("Samples with correlation between technical replicates <", samp.filt.thresh,sep="")
									write.table(bad_sample_names, file=feat.eval.outfile,sep="\t", row.names=FALSE)
								}
								
								bad_list={}
								if(length(bad_samples)>0)
								{
									for(n1 in 1:length(bad_samples))
									{	
										if(bad_samples[n1]>1)
										{
											bad_samples[n1]=bad_samples[n1]+(bad_samples[n1]-1)*(num_replicates-1)
										}
											
									}
									for(n1 in 1:num_replicates)
									{
										bad_list<-c(bad_list,(bad_samples+n1-1))
									}
									bad_list<-bad_list[order(bad_list)]
									
								}
								if(i>1){
										parent_bad_list<-intersect(parent_bad_list,bad_list)
									}
									else{
									    
									    parent_bad_list<-bad_list

									}
									
								
								#curdata<-cbind(curdata,feat.eval.result[which(as.numeric(feat.eval.result$median)<=feat.filt.thresh),])
								
								
								
							
                                  }
				  else
				  {
					  print("*********skipping feature evaluataion as only one replicate is present******")
				  }
                                
	       	}
	       	
		file_name="parameter_filenames.txt"
		pfile_name=paste(xMSanalyzer.outloc,"/",file_name,sep="")	
		write.table(pfilenames,file=pfile_name,sep="\t",row.names=FALSE)
		
		cat("\n")
		 print("*********Stage 2: Quality based filtering of results from each parameter setting based on sample and feature checks********")
		 cat("\n")
		  for(i in 1:length(alignmentresults))
                          {
                                  file_name=sapply(strsplit(alignmentresults[i],".txt"),head)
                                  feat.eval.file=paste(subdir2,"/",file_name,"cor",samp.filt.thresh,fileroot,feat.filt.thresh,deltamzminmax.tol,"ppmmzrangefiltereddata.txt",sep="")	
                                  #data_rpd[[i]]=read.table(feat.eval.file,header=TRUE)
				  
				  curdata<-data_rpd[[i]]
				 
				  feat.eval.result.mat<-feateval_list[[i]]
				  if(length(parent_bad_list)>0){
					  
					 
					  curdata<-curdata[,-c(parent_bad_list+4)]
					  #maxint<-apply(curdata[,-c(1:4,((dim(curdata)[2]-6):dim(curdata)[2]))],1,max)
					  maxint<-apply(curdata[,-c(1:4)],1,max)
					  badfeats<-which(maxint==0)
					  if(length(badfeats)>0){
						curdata<-curdata[-c(badfeats),]
						
						
						
						  feat.eval.result.mat<- feat.eval.result.mat[-c(badfeats),]
						  feateval_list[[i]]<-feat.eval.result.mat
					  }
					
					}
				 data_rpd[[i]]<-curdata[which(as.numeric(feat.eval.result.mat$median)<=feat.filt.thresh),]
				 
				
				 feat.eval.outfile=paste(subdir2,"/",file_name,"cor",samp.filt.thresh,fileroot,feat.filt.thresh,"filtereddata.txt",sep="")	
								
				 #write results
				 write.table(data_rpd[[i]], feat.eval.outfile,sep="\t", row.names=FALSE)
				 
				
			 }
                          ###########################################
                          #4) Merge two or more parameter settings
                          cat("\n")
                          print("*******Stage 3a: Merging features detected at different parameter settings**********")
							cat("\n")                         
                          num_pairs=1
                          finalres={}
                          rnames={}
                          
			  
			  for(i in 1:length(alignmentresults))
                          {
                                  file_name=sapply(strsplit(alignmentresults[i],".txt"),head)
                                  feat.eval.file=paste(subdir2,"/",file_name,"cor",samp.filt.thresh,fileroot,feat.filt.thresh,"filtereddata.txt",sep="")	
                                  #data_rpd[[i]]=read.table(feat.eval.file,header=TRUE)
				 

				  a1=sapply(strsplit(as.character(alignmentresults[i]), "\\_pres"), head, n=2)[2]
                                  minpres=sapply(strsplit(as.character(a1), "\\_run"), head, n=2)[1]
                                  minpres=as.numeric(minpres)
                                  a2=sapply(strsplit(as.character(alignmentresults[i]), "\\_run"), head, n=2)[2]
                                  minrun=sapply(strsplit(as.character(a2), "\\_"), head, n=2)[1]
                                  minrun=as.numeric(minrun)
                                  p1=paste(minrun,"_",minpres,sep="")
					    			

                                  for(j in i:length(alignmentresults))
                                  {
                                  	bool_num<-1
                                  	if(i==j){
                                  	if(length(alignmentresults)>1){
                                  		bool_num<-0
                                  	}
                                  	else{
                                  		bool_num<-1
                                  		}
                                  	}
                	#if(i!=j)
                	 if(bool_num==1)
				    {
		                          file_name=sapply(strsplit(alignmentresults[j],".txt"),head)
					  feat.eval.file=paste(subdir2,"/",file_name,"cor",samp.filt.thresh,fileroot,feat.filt.thresh,"filtereddata.txt",sep="")	
                                          #data_rpd[[j]]=read.table(feat.eval.file,header=TRUE)

                                          a1=sapply(strsplit(as.character(alignmentresults[j]), "\\_pres"), head, n=2)[2]
                                          minpres=sapply(strsplit(as.character(a1), "\\_run"), head, n=2)[1]
                                          minpres=as.numeric(minpres)
                                          a2=sapply(strsplit(as.character(alignmentresults[j]), "\\_run"), head, n=2)[2]
                                          minrun=sapply(strsplit(as.character(a2), "\\_"), head, n=2)[1]
                                          minrun=as.numeric(minrun)
                                          p2=paste(minrun,"_",minpres,sep="")
                                         # if(i!=j)
                                          {
                                                  p1_p2=paste("p",i,"_U_","p",j,sep="")
                                          }
                                          
                                          
					 feat.eval.A<-feateval_list[[i]]
					 feat.eval.B<-feateval_list[[j]]
					 feat.eval.A<-feat.eval.A[which(as.numeric(feat.eval.A$median)<=feat.filt.thresh),]
					 feat.eval.B<-feat.eval.B[which(as.numeric(feat.eval.B$median)<=feat.filt.thresh),]
					 
						print(paste("Number of good quality features from setting ",i,":", dim(data_rpd[[i]])[1],sep=": "))
						if(i!=j){
						print(paste("Number of good quality features from setting ",j,":",dim(data_rpd[[j]])[1],sep=": "))
						}
					 
                     data_m=merge.Results(data_rpd[[i]],data_rpd[[j]],feat.eval.A,feat.eval.B,max.mz.diff,max.rt.diff,merge.eval.pvalue,alignment.tool="apLCMS",
					 numnodes=numnodes,mult.test.cor,mergecorthresh,missingvalue)
					
					numcols<-dim(data_m)[2]

					 data_m_int<-data_m[,-c(1:5,(numcols-7):numcols)]
					 numsamps<-dim(data_m_int)[2]/num_replicates
					 
					 if(is.na(minexp.pct)==FALSE){
					minexp<-minexp.pct*dim(data_m_int)[2]
					
					if(length(which(data_m[,5]>=minexp))>0){
					data_m<-data_m[which(data_m[,5]>=minexp),]
					}else{
						stop(paste("No features have non-missing value in ",minexp, " samples",sep=""))
						}
					
					
					}
					
					
					 numcols<-dim(data_m)[2]

					 data_m_int<-data_m[,-c(1:5,(numcols-7):numcols)]
					 numsamps<-dim(data_m_int)[2]/num_replicates
					 
					 
					
					 
					maxint<-apply(data_m_int,1,function(x){max(x,na.rm=TRUE)})
					
					#numpeaks<-lapply(1:dim(data_m_int)[1],function(j){length(which(is.na(data_m_int[j,])==FALSE))})
				
					#if(length(which(data_m$numgoodsamples>minexp))>0){
						
					
					
					union_list[[num_pairs]]<-data_m[,c(1:5)]
					
					union_list[[num_pairs]]<-cbind(union_list[[num_pairs]],data_m$numgoodsamples)
					
					union_list[[num_pairs]]<-cbind(union_list[[num_pairs]],data_m$median)
					
					union_list[[num_pairs]]<-cbind(union_list[[num_pairs]],data_m$Qscore)
					
					
					union_list[[num_pairs]]<-cbind(union_list[[num_pairs]],maxint)
					
					union_list[[num_pairs]]<-cbind(union_list[[num_pairs]],data_m_int)
					
					featinfo<-colnames(data_m[,c(1:4)])
					cnames<-colnames(data_m_int)
					
					merge.res.colnames<-c(featinfo,"NumPres.All.Samples","NumPres.Biological.Samples",paste("median",fileroot,sep=""),"Qscore","Max.Intensity",cnames)
					colnames(union_list[[num_pairs]])<-as.character(merge.res.colnames)
					
					finalname=paste("apLCMSfeaturetable_at_", p1_p2,"cor",samp.filt.thresh,fileroot,feat.filt.thresh,".txt",sep="")
					  
                                          #Output merge results
                                          write.table(union_list[[num_pairs]],file=paste(subdir3,"/",finalname,sep=""), sep="\t",row.names=FALSE)
				
					
					  curres={}
					   curres=cbind(curres, p1_p2)
                                          curres=cbind(curres, dim(union_list[[num_pairs]])[1])
                                          curres=cbind(curres, mean(as.numeric(union_list[[num_pairs]][,7])))
                                          curres=cbind(curres, mean(as.numeric(data_m$Qscore)))
                                          
                                          curscore<-(dim(union_list[[num_pairs]])[1]-(scoreweight*mean(as.numeric(union_list[[num_pairs]][,7]))))
                                          
                                          if(curscore>bestscore){
                                          	
                                          	bestscore<-curscore
                                          	best_i<-num_pairs
                                          	best_pair<-p1_p2
                                          }
                                          
                                          
                                          curres=cbind(curres,curscore)
                                          
					  curres<-as.data.frame(curres)
                                          
                                          finalres=rbind(finalres,curres)
                                          

					 		 num_pairs=num_pairs+1
                       }
			      }				  
                          
                   }
                          finalres<-as.data.frame(finalres)
                         

			 
			  colnames(finalres)<-c("Parameter Combination", "Number of Features", "median PID/CV between sample replicates","mean Qscore (Quality score)", "Parameter score")
		      write.table(finalres,file=paste(xMSanalyzer.outloc,"/Stage3a/apLCMS_with_xMSanalyzer_merge_summary.txt",sep=""), sep="\t", row.names=FALSE)
		          
		          
		          


    
    print("Most optimal feature setting:")
    print(best_pair)
      cat("\n")
    #########################################################################
    
    
    cat("\n")
      print(paste("********Stage 3b: Generating final (pre-batcheffect correction) untargeted and targeted feature tables using ",best_pair," results******",sep=""))
      cat("\n")
      #pdf(file=paste("subdir4a/Stage4a_QC_plots.pdf",sep=""))
      
      pdf(file=paste(subdir4a,"/Stage4a_QC_plots.pdf",sep=""))
      
 	d1<-union_list[[best_i]]
 	

  
   
   			
			
				max_numzeros<-dim(d1)[2]*1

	time_thresh<-NA
	
	if(is.na(void.vol.timethresh)==FALSE){
		dfirst15<-d1[which(d1[,2]<void.vol.timethresh),]
	if(length(dfirst15)>0){
		ind1<-which(dfirst15[,9]==max(dfirst15[,9]))
		
		time_thresh<-dfirst15[ind1,2]
		
		time_thresh<-time_thresh-(0.30*time_thresh)
		
		dfirst15<-dfirst15[order(dfirst15[,2]),]

		
		plot(dfirst15[,2],dfirst15[,9],col="red",xlab="Time (s)", ylab="Max intensity across all samples", main=paste("Estimated void volume time: ",time_thresh," s \n (using all data)",sep=""),cex.main=0.7)
		abline(v=time_thresh,col=4,lty=3)
			
		d1<-d1[which(d1$time>=time_thresh),]
		    
		 print("Estimated void volume time cutoff")
		 print(time_thresh)
		    
		    }
		d1_int<-round(d1[,-c(1:9)],0)
		d1<-cbind(d1[,c(1:9)],d1_int)
		rm(d1_int)
			
		
		#sfname<-paste(xMSanalyzer.outloc,"/stage5/feature_table_",best_pair,"_cor",samp.filt.thresh,fileroot,feat.filt.thresh,"filtered.txt",sep="")
		finalname<-paste("featuretable_",best_pair,"_cor",samp.filt.thresh,fileroot,feat.filt.thresh,"_voidtimefilt.txt",sep="")
		
		write.table(d1,file=paste(subdir3b,"/",finalname,sep=""),sep="\t",row.names=FALSE)
	   

	}else{
		
		
		finalname<-paste("featuretable_",best_pair,"_cor",samp.filt.thresh,fileroot,feat.filt.thresh,".txt",sep="")
		
		write.table(d1,file=paste(subdir3b,"/",finalname,sep=""),sep="\t",row.names=FALSE)
	   

		
		}

 cat("\n")
   print(paste("********Stage 4a: Performing QC checks using ",best_pair," results******",sep=""))
         
      cat("\n")
      Sys.sleep(1)
	
    
		
		print("Dim data after void time filtering")
		print(dim(d1))

hist(d1$NumPres.All.Samples,main="Histogram NumPres.All.Samples \n (using all data)",col="orange",xlab="Number of samples (including replicates) \n with non-zero intensity values", ylab="Number of features",cex.main=0.7)
			
			hist(d1$NumPres.Biological.Samples,main="Histogram NumPres.Biological.Samples \n (using all data)",col="orange",xlab="Number of biological samples \n with non-zero intensity values", ylab="Number of features",cex.main=0.7)
			#h1<-hist(d1$median_CV,seq(0,max(d1$median_CV,na.rm=TRUE),10),main="Histogram median CV \n (using all data)",col="orange",xlab="median CV%", ylab="Number of features",cex.main=0.7)
			
			#h1<-hist(d1$median_CV,breaks=seq(0,max(d1$median_CV,na.rm=TRUE)+10,10),main="Histogram median CV \n (using all data)",col="orange",xlab="median CV%", ylab="Number of features",cex.main=0.7)
			
			if(num_replicates>2){
			h1<-hist(d1$median_CV,breaks=seq(0,max(d1$median_CV,na.rm=TRUE)+10,10),main="Histogram median CV \n (using all data)",col="orange",xlab="median CV%", ylab="Number of features",cex.main=0.7)
		}else{
			h1<-hist(d1$median_PID,breaks=seq(0,max(d1$median_PID,na.rm=TRUE)+10,10),main="Histogram median CV \n (using all data)",col="orange",xlab="median CV%", ylab="Number of features",cex.main=0.7)
			}	
			
				
			hist(d1$Qscore,main="Histogram Qscore \n (using all data)",col="orange",xlab="Quality score", ylab="Number of features",cex.main=0.7)
			
			lab_text<-paste(h1$breaks,"-",h1$breaks+10,sep="")
			#pie(h1$counts,col=rainbow(length(h1$counts)),labels=lab_text,main="Pie chart of median CVs (%) \n using all features",cex.main=0.7)
			
			if(num_replicates>2){
			pie(h1$counts,col=rainbow(length(h1$counts)),labels=lab_text,main=paste("Pie chart of median CVs (%) \n using all features\n; average=",mean(d1$median_CV),sep=""),cex.main=0.7)
			}else{
				pie(h1$counts,col=rainbow(length(h1$counts)),labels=lab_text,main=paste("Pie chart of median PIDs (%) \n using all features\n; average=",mean(d1$median_PID),sep=""),cex.main=0.7)
			
				}
			d2<-d1[order(d1$time),]
			d2<-as.data.frame(d2)
			#d2<-apply(d2,1,as.numeric)
			
			plot(d2$time,d2$Max.Intensity,main="TIC \n (using all data)",col="orange",xlab="Time (s)",ylab="Max intensity across all samples/profiles",cex.main=0.7)
			
			plot(d2$time,d2$mz,main="m/z vs Time \n (using all data)",col="orange",xlab="Time (s)",ylab="m/z",cex.main=0.7)
			
			plot(d2$time,d2$Qscore,main="Qscore vs Time \n (using all data)",col="orange",xlab="Time (s)",ylab="Qscore",cex.main=0.7)
			
			plot(d2$mz,d2$Qscore,main="Qscore vs m/z \n (using all data)",col="orange",xlab="m/z",ylab="Qscore",cex.main=0.7)
			
		
			
	
	      
	data_m<-d1[,-c(1:9)]
   
    if(replacezeroswithNA==TRUE){
	data_m<-replace(as.matrix(data_m),which(data_m==0),NA)
	}
	
	
	counna<-apply(data_m,1,function(x){length(which(is.na(x)==TRUE))})

	maxzeros<-1*dim(data_m)[2]
	
	if(length(which(counna<maxzeros))){
		data_m<-data_m[which(counna<maxzeros),]
	}

    maxint<-apply(data_m,1,function(x){max(x,na.rm=TRUE)})
    
   	maxint_ord<-order(maxint,decreasing=TRUE)
   	#[maxint_ord[1:5000]
    X<-t(data_m) #[maxint_ord[1:2000],])
    X<-replace(as.matrix(X),which(is.na(X)==TRUE),0)
    
    tic.eval(d1[,-c(1:9)],outloc=subdir4a)
    
    feat.eval.result=evaluate.Features(d1[,-c(5:9)], numreplicates=num_replicates,min.samp.percent=0.6,alignment.tool="apLCMS",impute.bool=TRUE)
										cnames=colnames(feat.eval.result)
		
				Sys.sleep(1)
			


    
    Sys.sleep(1)
    
s1<-"Stage 1 results: QC evaluation of invidual parameters from apLCMS"
s2<-"Stage 2 results: filtered results from each paramter setting based on sample and feature quality (CV within replicates) checks"
s3<-"Stage 3a results: merged results using stage 2 filtered files"
s3b<-"Stage 3b results: Final untargeted and targeted feature tables"
s4a<-"Stage 4a results: QC evaluation of targeted and untargeted data before batch-effect correction"

sm<-rbind(s1,s2,s3,s3b,s4a)

targeted_feat_raw<-eval.target.mz(dataA=d1[,-c(3:9)],refMZ=stddata,feature.eval.result=feat.eval.result,mzthresh=refMZ.mz.diff,timethresh=refMZ.time.diff,outloc=subdir4a,folderheader="raw")
		

    if(is.na(sample_info_file)==FALSE)
    {
    	#sampleid_mapping<-read.table(sample_info_file,sep="\t",header=TRUE)
    	
    	sampleid_mapping[,1]<-gsub(sampleid_mapping[,1],pattern=filepattern,replacement="")
		cnames<-colnames(data_m)
		
		cnames<-gsub(cnames,pattern=filepattern,replacement="")
		colnames(data_m)<-cnames


    	batch_inf<-sampleid_mapping[,3]
    	
    	batch_labels<-as.factor(sampleid_mapping[,3])
    
    	l1<-levels(batch_labels)

		    pca.eval(X=X,samplelabels=batch_labels,filename="raw",ncomp=5,center=TRUE,scale=TRUE,legendlocation="topright",legendcex=0.5,outloc=subdir4a)
			 
			Sys.sleep(1)
			
			
			
			
			   		    
    		if(dim(sampleid_mapping)[2]<4){
   		    mod<-rep(1,dim(data_m)[2])
   		    }else{
   		    		mod<-sampleid_mapping[,-c(1:3)]
   		    	}
    

  	
    		dev.off()
    		dev.off()
    		
    		
    		Sys.sleep(5)
    
		     #######################################################################
		    cat("\n")
		     print(paste("Stage 4b: ComBat processing and post-ComBat QC evaluation using ",best_pair," results",sep=""))
		    cat("\n")
		    #pdf("Stage4b_QC_plots.pdf")
		   # pdf(file=paste("subdir4b/Stage4b_QC_plots.pdf",sep=""))
		    pdf(file=paste(subdir4b,"/Stage4b_QC_plots.pdf",sep=""))
		    ##################################################################
		
		    adjdata<-try(sva::ComBat(dat=data_m,batch=batch_inf,mod=mod,par.prior=TRUE),silent=TRUE)
		    
		    if(is(adjdata,"try-error")){
		 
		
				data_m1<-cbind(d1[,c(1:4)],data_m)
				
				#print(adjdata)
				
				print("sva::ComBat caused an error.")
				print(adjdata)
				
				print("Too many missing values can cause this error. Trying MetabComBat...")
				adjdata<-MetabComBat(dat=data_m1,saminfo=sampleid_mapping,par.prior=T,filter=F,write=F,prior.plots=F)
		    	adjdata<-adjdata[,-c(1:4)]
		    }
    
    				  
		    
		    maxint<-apply(adjdata,1,function(x){max(x,na.rm=TRUE)})
		    
		    adjdata2<-replace(as.matrix(adjdata),which(is.na(adjdata)==TRUE),0)
		    adjdata2<-replace(as.matrix(adjdata2),which(adjdata2<0),0)
		    
		    adjdata2<-cbind(d1[,c(1:4)],adjdata2)
		    
		    
		    
		    expression_xls<-paste(subdir4b,"/ComBat_corrected_",best_pair,"_feature_table.txt",sep="")
		    
		    
		   							
			cnames=colnames(feat.eval.result)
										
			feat.eval.result<-apply(feat.eval.result,2,as.numeric)
			feat.eval.result<-as.data.frame(feat.eval.result)
			feat.eval.result.mat=cbind(adjdata2[,c(1:4)],feat.eval.result)  
			
			numpeaks<-apply(adjdata2[,-c(1:4)],1,countpeaks)						
		    
		    maxint<-apply(adjdata2[,-c(1:4)],1,function(x){max(x,na.rm=TRUE)})
		    
		    adjdata2<-cbind(adjdata2[,c(1:4)],numpeaks,feat.eval.result.mat$numgoodsamples,feat.eval.result.mat$median,feat.eval.result.mat$Qscore,maxint,
		    adjdata2[,-c(1:4)])
		    
		    colnames(adjdata2)<-colnames(d1)
		    setwd(subdir4b)
		    write.table(adjdata2,file=expression_xls,sep="\t",row.names=FALSE)
		    
		   # X<-t(adjdata2[maxint_ord[1:2000],-c(1:9)])
		   
		    feat.eval.result=evaluate.Features(adjdata2, numreplicates=num_replicates,min.samp.percent=0.6,alignment.tool="apLCMS",impute.bool=TRUE)
			
		    
		     tic.eval(adjdata2[,-c(1:9)],outloc=subdir4b)
    
		   
		     X<-t(adjdata2[,-c(1:9)])
		   
		     
		    X<-as.matrix(X)
		    X<-replace(as.matrix(X),which(is.na(X)==TRUE),0)
		    pca.eval(X,samplelabels=batch_labels,filename="ComBat",ncomp=5,center=TRUE,scale=TRUE,legendlocation="topright",legendcex=0.5,outloc=subdir4b)
			    
		    
		    #eval.reference.metabs(dataA=adjdata2,stdData=stddata,mzthresh=10,outloc=getwd(),folderheader="ComBat")
		    
			targeted_feat_combat<-eval.target.mz(dataA=adjdata2[,-c(3:9)],refMZ=stddata,feature.eval.result=feat.eval.result,mzthresh=refMZ.mz.diff,timethresh=refMZ.time.diff,outloc=subdir4b,folderheader="ComBat")
			
		s4b<-"Stage 4b results: Batch-effect evaluation post ComBat and QC evaluation of targeted data post batch-effect correction"
		sm<-rbind(sm,s4b)
	
		}else{
			
			adjdata2=NULL
			targeted_feat_combat<-NULL
			}
    
			dev.off()
			
					 metlin.res={}
                      kegg.res={}
                      
					 #length(union_list[[num_pairs]]$mz
					 if(is.na(adduct.list)==FALSE){
					
					cat("\n")
					  print("*********Stage 5: Mapping m/z values to known metabolites using KEGG*********")
					 cat("\n")
					 
                     	annot.res<-feat.batch.annotation.KEGG(adjdata2,mz.tolerance.dbmatch,adduct.list,subdir5, numnodes=numnodes,syssleep=syssleep)
					    annotres_list<-annot.res
					    s5<-"Stage 5 results: Annotation of features"
					    sm<-rbind(sm,s5)
                     }
				

		          
               		  print("*************Processing complete**********") 

                }
                else
                {
                        stop(paste("No files exist in",apLCMS.outloc, "Please check the input value for cdfloc", sep=""))
                }
                
                suppressWarnings(sink(file=NULL))

sm<-as.data.frame(sm)
colnames(sm)<-"Output_description"
setwd(xMSanalyzer.outloc)
write.table(sm,file="Readme.txt",sep="\t",row.names=FALSE)
                                print(paste("Processing is complete. Program results can be found at: ",xMSanalyzer.outloc,sep=""))
		
        #return(list("mergeresult"=union_list, "apLCMSres"=data_rpd_all,"apLCMSresfilt"=data_rpd))
	return(list("apLCMS.merged.res"=union_list, "apLCMS.ind.res"=data_rpd_all,"apLCMS.ind.res.filtered"=data_rpd,"final.feat.table.annot"=annotres_list, "feat.eval.ind"=feateval_list, "sample.eval.ind"=rsqres_list,"final.feat.table.raw"=d1,"final.feat.table.combat"=adjdata2,
	"final.targeted.feat.table.raw"=targeted_feat_raw,"final.targeted.feat.table.combat"=targeted_feat_combat))
	
	
}
