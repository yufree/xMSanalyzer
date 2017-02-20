countpeaks <-
function(intensity_vec, missing_val=0){
                
						if(is.na(missing_val)==TRUE){
							zeroind<-which(is.na(intensity_vec)==FALSE)
						}else{
							if(missing_val==0){
								zeroind<-which(intensity_vec>0)
							}
							else{
							 stop(paste("Invalid value for \"missing_val\". Please use either \"NA\" or \"0\"", sep=""))
							}
						}
						

						return(length(zeroind))
					}
