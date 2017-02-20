get_deltappm <-
function(x){
	if(is.na(x[1])==FALSE){
	return(10^6*(abs(x[1]-x[2])/x[1]))
	}else{
	return(0)
	}
}
