Hothorn_gnrt <- function(N=200, Model.type = 2, Dist = "Exponential", censor.type=1, truncation = 2){
	
	Data <- as.data.frame(matrix(NA,N,11))
	names(Data)<-c("X1","X2","X3","X4","X5","X6","Start","T","C","Obs","Event")
	Count = 0
	
	while(Count < N){
		x1 <- runif(1,0,1)
		x2 <- sample(c(0,1),1)
		x3 <- sample(c(0,1),1)
		x4 <- runif(1,0,1)
		x5 <- runif(1,0,1)
		x6 <- sample(c(0,1),1)
		L  <- runif(1,0,truncation)

		if(Model.type == 3){
			Param <- -(cos((x1+x2)*pi)+sqrt(x1+x2))
		}else if(Model.type == 2){
			Param <- -x1-x2
		}else {
			stop("Wrong model type: It's either 2(second setup) or 3 (third setup)")
		}

		if(Dist == "Exponential"){
			t <- rexp(1,exp(Param))
		}else if(Dist == "Weibull-I"){
			t <- rweibull(1, shape = 2, scale = 10*exp(Param))
		}else if(Dist == "Weibull-D"){
			t <- rweibull(1, shape = 0.5, scale = 5*exp(Param))
		}else{
			print("Error: Wrong distribution!");
			return(0);
		}

		if( t > L){
			Count = Count + 1

			Data[Count,"X1"] <- x1
			Data[Count,"X2"] <- x2
			Data[Count,"X3"] <- x3
			Data[Count,"X4"] <- x4
			Data[Count,"X5"] <- x5
			Data[Count,"X6"] <- x6
			Data[Count,"Start"] <- L
			Data[Count,"T"] <- t
			
			if(Dist == "Exponential"){
				if(censor.type==1){
					Data[Count,"C"] <- L + rexp(1,1/13)
				}else if(censor.type==2){
					Data[Count,"C"] <- L + rexp(1,1/3)
				}
			}else if(Dist == "Weibull-D"){
				 if(censor.type==1){
					Data[Count,"C"] <- L + rexp(1,1/23)
				}else if(censor.type==2){
					Data[Count,"C"] <- L + rexp(1,1/4.2)
				}

			}else if(Dist == "Weibull-I"){
				if(censor.type==1){
					Data[Count,"C"] <- L + rexp(1,1/14)
				}else if(censor.type==2){
					Data[Count,"C"] <- L + rexp(1,1/3.5)
				}
			}

			if( t <= Data[Count,"C"]){
				Data[Count,"Obs"] <- t
				Data[Count,"Event"] <- 1
			}else{
				Data[Count,"Obs"] <- Data[Count,"C"]
				Data[Count,"Event"] <- 0
			}

		}#end if
	}#end while 

	return(Data)
}

###############################

Hothorn_gnrt_test <- function(N=200, Model.type = 2, Dist = "Exponential"){
	Data <- as.data.frame(matrix(NA,N,8))
	names(Data)<-c("X1","X2","X3","X4","X5","X6","Obs","Event")
	
	Data$X1 <- runif(N,0,1)
	Data$X2 <- sample(c(1,0),N,replace=TRUE)
	Data$X3 <- sample(c(1,0),N,replace=TRUE)
	Data$X4 <- runif(N,0,1)
	Data$X5 <- runif(N,0,1)
	Data$X6 <- sample(c(1,0),N,replace=TRUE)

	for(k in 1:N){
		if(Model.type == 3){
			Param <- -(cos((Data[k,"X1"]+Data[k,"X2"])*pi)+sqrt(Data[k,"X1"]+Data[k,"X2"]))
		}else if(Model.type == 2){
			Param <- -Data[k,"X1"]-Data[k,"X2"]
		}else {
			stop("Wrong model type: It's either 2(second setup) or 3 (third setup)")
		}

		if(Dist == "Exponential"){
			Data[k,"Obs"] <- rexp(1,exp(Param))
		}else if(Dist == "Weibull-I"){
			Data[k,"Obs"] <- rweibull(1, shape = 2, scale = 10*exp(Param))
		}else if(Dist == "Weibull-D"){
			Data[k,"Obs"] <- rweibull(1, shape = 0.5, scale = 5*exp(Param))
		}else{
			print("Error: Wrong distribution!");
			return(0);
		}
	}##end of for loop
	Data$Event <- 1
	return(Data)
}
	



