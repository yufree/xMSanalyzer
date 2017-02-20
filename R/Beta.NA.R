Beta.NA <-
function(y,X){
	des=X[!is.na(y),]
	#print(dim(des))
	sum_des<-apply(des,2,sum)
	
	#write.table(des,file="current_des1.txt",sep="\t",row.names=TRUE)
	
	bad_batches<-which(sum_des==0)
	if(length(bad_batches)>0){
	des<-des[,-bad_batches]
	}
	#write.table(des,file="current_des2.txt",sep="\t",row.names=TRUE)
	y1=y[!is.na(y)]
	B <- solve(t(des)%*%des)%*%t(des)%*%y1
	#print(dim(B))
	
	Ball<-matrix(0,nrow=length(sum_des),ncol=1)
	if(length(bad_batches)>0){
	Ball[-bad_batches,]<-B
	}else{
	Ball<-B
	}
	#print(Ball[1:10,])
	Ball
	}
