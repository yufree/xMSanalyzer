MetabComBat <-
function(dat, saminfo, type='txt', write=T, covariates='all', par.prior=T, filter=F, skip=0, prior.plots=T){
	#debug: expression_xls='exp.txt'; sample_info_file='sam.txt'; type='txt'; write=T; covariates='all'; par.prior=T; filter=F; skip=0; prior.plots=T
	#cat('Reading Sample Information File\n')
	#saminfo <- read.table(sample_info_file, header=T, sep='\t',comment.char='')
	if(sum(colnames(saminfo)=="Batch")!=1){stop('ERROR: Sample Information File does not have a Batch column!');
		return('ERROR: Sample Information File does not have a Batch column!')}
		
		if (skip>0){
              geneinfo <- as.matrix(dat[,1:skip])
              dat <- dat[,-c(1:skip)]
	}else{geneinfo=dat[,c(1:4)]; dat<-dat[,-c(1:4)]}
        #print(geneinfo[1:4])
        print(dat[1:2,1:3])
	
	if(filter){
		nfeatures <- nrow(dat)
		col <- ncol(dat)/2
		present <- apply(dat, 1, filter.absent, filter)
		dat <- dat[present, -(2*(1:col))]
		if (skip>0){geneinfo <- geneinfo[present,]}
		cat('Filtered features absent in more than',filter,'of samples. features remaining:',nrow(dat),'; features filtered:',nfeatures-nrow(dat),'\n')
		}

	if(any(apply(dat,2,mode)!='numeric')){stop('ERROR: Array expression columns contain non-numeric values! (Check your .xls file for non-numeric values and if this is not the problem, make a .csv file and use the type=csv option');return('ERROR: Array expression columns contain non-numeric values! (Check your .xls file for non-numeric values and if this is not the problem, make a .csv file and use the type=csv option)')}
	
	#cnames1<-colnames(dat)
	tmp <- match(colnames(dat),saminfo[,1])
	#tmp<-which(cnames1%in%saminfo[,1])
	#tmplen<-length(tmp)
	#if(tmplen!=length(cnames1)){return('ERROR: Sample Information File and Data Array Names are not the same!')}
	if(any(is.na(tmp))){stop('ERROR: Sample Information File and Data Array Names are not the same!');return('ERROR: Sample Information File and Data Array Names are not the same!')}
	#tmp1 <- match(saminfo[,1],colnames(dat))
	#saminfo <- saminfo[tmp1[!is.na(tmp1)],]
	saminfo <- saminfo[tmp,]  ## Bug fixed 01/04/2011		

	if(any(covariates != 'all')){saminfo <- saminfo[,c(1:2,covariates)]}
	design <- design.mat(saminfo)	


	batches <- list.batch(saminfo)
	n.batch <- length(batches)
	n.batches <- sapply(batches, length)
	n.array <- sum(n.batches)
	
	## Check for missing values
	NAs = any(is.na(dat))
	if(NAs){cat(c('Found',sum(is.na(dat)),'Missing Data Values\n'),sep=' ')}
        #print(dat[1:2,])
	##Standardize Data across features
	cat('Standardizing Data across features\n')
	if (!NAs){B.hat <- solve(t(design)%*%design)%*%t(design)%*%t(as.matrix(dat))
	
	}else{B.hat=apply(dat,1,Beta.NA,design)
	#print(length(B.hat))
	
	#print(dim(B.hat[[1]]))
	
	#save(B.hat,file="bhat.Rda")
	B.hat_orig<-B.hat
	
	
	numfeats<-length(B.hat)/length(n.batches)
	#print(numfeats)
	#print(n.batches)
	dim(B.hat)<-dim(matrix(0,nrow=length(n.batches),ncol=numfeats))
	
	} #Standarization Model
	
	#rows: number of batches; columns: number of metabs
	print("standardization done")
	print(dim(B.hat))
	
	
	grand.mean <- t(n.batches/n.array)%*%B.hat[1:n.batch,]
	if (!NAs){var.pooled <- ((dat-t(design%*%B.hat))^2)%*%rep(1/n.array,n.array)}else{var.pooled <- apply(dat-t(design%*%B.hat),1,var,na.rm=T)}

	stand.mean <- t(grand.mean)%*%t(rep(1,n.array))
	if(!is.null(design)){tmp <- design;tmp[,c(1:n.batch)] <- 0;stand.mean <- stand.mean+t(tmp%*%B.hat)}	
	s.data <- (dat-stand.mean)/(sqrt(var.pooled)%*%t(rep(1,n.array)))

	##Get regression batch effect parameters
	cat("Fitting L/S model and finding priors\n")
	batch.design <- design[,1:n.batch]
	if (!NAs){gamma.hat <- solve(t(batch.design)%*%batch.design)%*%t(batch.design)%*%t(as.matrix(s.data))}else{gamma.hat=apply(s.data,1,Beta.NA,batch.design)}
	delta.hat <- NULL
	for (i in batches){
		delta.hat <- rbind(delta.hat,apply(s.data[,i], 1, var,na.rm=T))
		}

	delta.hat<-replace(delta.hat,which(is.na(delta.hat)==TRUE),1)
	
			
	##Find Priors
	gamma.bar <- apply(gamma.hat, 1, mean)
	t2 <- apply(gamma.hat, 1, var)
	a.prior <- apply(delta.hat, 1, aprior)
	b.prior <- apply(delta.hat, 1, bprior)

	
	##Plot empirical and parametric priors

	if (prior.plots & par.prior){
		par(mfrow=c(2,2))
		tmp <- density(gamma.hat[1,])
		plot(tmp,  type='l', main="Density Plot")
		xx <- seq(min(tmp$x), max(tmp$x), length=100)
		lines(xx,dnorm(xx,gamma.bar[1],sqrt(t2[1])), col=2)
		qqnorm(gamma.hat[1,])	
		qqline(gamma.hat[1,], col=2)	
	
		tmp <- density(delta.hat[1,])
		invgam <- 1/rgamma(ncol(delta.hat),a.prior[1],b.prior[1])
		tmp1 <- density(invgam)
		plot(tmp,  typ='l', main="Density Plot", ylim=c(0,max(tmp$y,tmp1$y)))
		lines(tmp1, col=2)
		qqplot(delta.hat[1,], invgam, xlab="Sample Quantiles", ylab='Theoretical Quantiles')	
		lines(c(0,max(invgam)),c(0,max(invgam)),col=2)	
		title('Q-Q Plot')
	}
	
	##Find EB batch adjustments

	gamma.star <- delta.star <- NULL
	if(par.prior){
		cat("Finding parametric adjustments\n")
		for (i in 1:n.batch){
			#delta.hat[i,]<-replace(delta.hat[i,],which(is.na(delta.hat[i,])==TRUE),1)
			
			

			temp <- it.sol(s.data[,batches[[i]]],gamma.hat[i,],delta.hat[i,],gamma.bar[i],t2[i],a.prior[i],b.prior[i])
			gamma.star <- rbind(gamma.star,temp[1,])
			delta.star <- rbind(delta.star,temp[2,])
			}
	}else{
		cat("Finding nonparametric adjustments\n")
		for (i in 1:n.batch){
			bad_cols<-which(is.na(delta.hat[i,])==TRUE)
			#delta.hat[i,]<-replace(delta.hat[i,],which(is.na(delta.hat[i,])==TRUE),1)
			
			
			temp <- int.eprior(as.matrix(s.data[,batches[[i]]]),gamma.hat[i,],delta.hat[i,])
			gamma.star <- rbind(gamma.star,temp[1,])
			delta.star <- rbind(delta.star,temp[2,])
			}
		}


	### Normalize the Data ###
	cat("Adjusting the Data\n")

	bayesdata <- s.data
	j <- 1
	for (i in batches){
		bayesdata[,i] <- (bayesdata[,i]-t(batch.design[i,]%*%gamma.star))/(sqrt(delta.star[j,])%*%t(rep(1,n.batches[j])))
		j <- j+1
		}

	bayesdata <- (bayesdata*(sqrt(var.pooled)%*%t(rep(1,n.array))))+stand.mean
	if(write){
		output_file <- paste('Adjusted_data_parprior',par.prior,'.xls',sep='_')
                 #print(geneinfo[1:2])
                 #print(bayesdata[1:2,1:4])
		 #cat(c(colnames(geneinfo),colnames(dat),'\n'),file=output_file,sep='\t')
		#suppressWarnings(write.table(cbind(geneinfo,formatC(as.matrix(bayesdata), format = "f")), file=output_file, sep="\t", quote=F,row.names=F,col.names=F,append=T))
                outdata <- cbind(ProbeID=geneinfo, bayesdata); write.table(outdata, file=output_file, sep="\t")
		cat("Adjusted data saved in file:",output_file,"\n")
		}else{return(cbind(geneinfo,bayesdata))}
	}
