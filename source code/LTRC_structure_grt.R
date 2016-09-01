source("bathtub.R")

LTRC.generate <- function(n=200, Dist = "Exponential", censor.type=1, truncation = 2){###censor.sup = "Wide"###
	
	Data <- as.data.frame(matrix(NA,n,8))
	names(Data)<-c("X1","X2","X3","Start","T","C","Obs","Event")
	Count = 0

	while (Count < n){
		x1 <- sample(1:5,1)
		x2 <- sample(c(1,2),1) 
		x3 <- runif(1,0,2)
		Ll = runif(1,0,truncation)


		if(Dist == "Bathtub"){
			if(x1 < 2.5){
				if(x2 == 1){
					t <- Bathtub(n=1, a = 0.01)
				}else{
					t <- Bathtub(n=1, a = 0.05)
				}
			}else{
				if(x3 <= 1){
					t <- Bathtub(n=1, a = 0.1)
				}else{
					t <- Bathtub(n=1, a = 0.7)
				}
			}
		}else if(Dist == "Exponential"){
			if(x1 < 2.5){
				if(x2 == 1){
					t <- rexp(1,0.1)
				}else{
					t <- rexp(1,0.23)
				}
			}else{
				if(x3 <= 1){
					t <- rexp(1,0.4)
				}else{
					t <- rexp(1,0.9)
				}
			}
		}else if(Dist == "Weibull-D"){
			if(x1 < 2.5){
				if(x2 == 1){
					t <- rweibull(1, 0.9, 7)
				}else{
					t <- rweibull(1, 0.9, 3)
				}
			}else{
				if(x3 <= 1){
					t <- rweibull(1, 0.9, 2.5)
				}else{
					t <- rweibull(1, 0.9, 1)
				}
			}
		}else if(Dist == "Weibull-I"){
			if(x1 < 2.5){
				if(x2 == 1){
					t <- rweibull(1, 3, 10)
				}else{
					t <- rweibull(1, 3, 6.2)
				}
			}else{
				if(x3 <= 1){
					t <- rweibull(1, 3, 4.3)
				}else{
					t <- rweibull(1, 3, 2)
				}
			}
		}else if(Dist == "Lognormal"){
			if(x1 < 2.5){
				if(x2 == 1){
					t <- rlnorm(1, meanlog = 2.0, sdlog = 0.3)
				}else{
					t <- rlnorm(1, meanlog = 1.7, sdlog = 0.2)
				}
			}else{
				if(x3 <= 1){
					t <- rlnorm(1, meanlog = 1.3, sdlog = 0.3)
				}else{
					t <- rlnorm(1, meanlog = 0.5, sdlog = 0.5)
				}
			}
		}

		if( t >= Ll){
			Count = Count + 1

			Data[Count,"X1"] <- x1
			Data[Count,"X2"] <- x2
			Data[Count,"X3"] <- x3
			Data[Count,"Start"] <- Ll
			Data[Count,"T"] <- t
			if(censor.type==1){
				if(Dist == "Bathtub"){
					Data[Count,"C"] <- Ll + rexp(1,1/14)
				}else if(Dist == "Exponential"){
					Data[Count,"C"] <- Ll + rexp(1,1/17)
				}else if(Dist == "Weibull-D"){
					Data[Count,"C"] <- Ll + rexp(1,1/15)
				}else if(Dist == "Weibull-I"){
					Data[Count,"C"] <- Ll + rexp(1,1/17)
				}else if(Dist == "Lognormal"){
					Data[Count,"C"] <- Ll + rexp(1,1/16)
				}
			}else{
				if(Dist == "Bathtub"){
					Data[Count,"C"] <- Ll + rexp(1,1/3.7)
				}else if(Dist == "Exponential"){
					Data[Count,"C"] <- Ll + rexp(1,1/4)
				}else if(Dist == "Weibull-D"){
					Data[Count,"C"] <- Ll + rexp(1,1/3.5)
				}else if(Dist == "Weibull-I"){
					Data[Count,"C"] <- Ll + rexp(1,1/4.5)
				}else if(Dist == "Lognormal"){
					Data[Count,"C"] <- Ll + rexp(1,1/4.4)
				}
			}
			
			if( t <= Data[Count,"C"]){
				Data[Count,"Obs"] <- t
				Data[Count,"Event"] <- 1
			}else{
				Data[Count,"Obs"] <- Data[Count,"C"]
				Data[Count,"Event"] <- 0
			}
		}
	}##end of while 

	Data$X4 <- sample(1:5,n,replace = TRUE)
	Data$X5 <- sample(c(1,2),n,replace = TRUE)
	Data$X6 <- runif(n,0,2)
	return(Data)
}

####==============================================================================

LTRC.generate.test <- function(n=200, Dist = "Exponential"){###censor.sup = "Wide"###
	
	Data <- as.data.frame(matrix(NA,n,5))
	names(Data)<-c("X1","X2","X3","Obs","Event")

	Data$X1 <- sample(1:5,n,replace=TRUE)
	Data$X2 <- sample(c(1,2),n,replace=TRUE)
	Data$X3 <- runif(n,0,2)

	for(k in 1:n){
		if(Dist == "Bathtub"){
			if(Data[k,"X1"] < 2.5){
				if(Data[k,"X2"] == 1){
					Data[k,"Obs"] <- Bathtub(n=1, a = 0.01)
				}else{
					Data[k,"Obs"] <- Bathtub(n=1, a = 0.05)
				}
			}else{
				if(Data[k,"X3"] <= 1){
					Data[k,"Obs"] <- Bathtub(n=1, a = 0.1)
				}else{
					Data[k,"Obs"] <- Bathtub(n=1, a = 0.7)
				}
			}
		}else if(Dist == "Exponential"){
			if(Data[k,"X1"] < 2.5){
				if(Data[k,"X2"] == 1){
					Data[k,"Obs"] <- rexp(1,0.1)
				}else{
					Data[k,"Obs"] <- rexp(1,0.23)
				}
			}else{
				if(Data[k,"X3"] <= 1){
					Data[k,"Obs"] <- rexp(1,0.4)
				}else{
					Data[k,"Obs"] <- rexp(1,0.9)
				}
			}
		}else if(Dist == "Weibull-D"){
			if(Data[k,"X1"] < 2.5){
				if(Data[k,"X2"] == 1){
					Data[k,"Obs"] <- rweibull(1, 0.9, 7)
				}else{
					Data[k,"Obs"] <- rweibull(1, 0.9, 3)
				}
			}else{
				if(Data[k,"X3"] <= 1){
					Data[k,"Obs"] <- rweibull(1, 0.9, 2.5)
				}else{
					Data[k,"Obs"] <- rweibull(1, 0.9, 1)
				}
			}
		}else if(Dist == "Weibull-I"){
			if(Data[k,"X1"] < 2.5){
				if(Data[k,"X2"] == 1){
					Data[k,"Obs"] <- rweibull(1, 3, 10)
				}else{
					Data[k,"Obs"] <- rweibull(1, 3, 6.2)
				}
			}else{
				if(Data[k,"X3"] <= 1){
					Data[k,"Obs"] <- rweibull(1, 3, 4.3)
				}else{
					Data[k,"Obs"] <- rweibull(1, 3, 2)
				}
			}
		}else if(Dist == "Lognormal"){
			if(Data[k,"X1"] < 2.5){
				if(Data[k,"X2"] == 1){
					Data[k,"Obs"] <- rlnorm(1, meanlog = 2.0, sdlog = 0.3)
				}else{
					Data[k,"Obs"] <- rlnorm(1, meanlog = 1.7, sdlog = 0.2)
				}
			}else{
				if(Data[k,"X3"] <= 1){
					Data[k,"Obs"] <- rlnorm(1, meanlog = 1.3, sdlog = 0.3)
				}else{
					Data[k,"Obs"] <- rlnorm(1, meanlog = 0.5, sdlog = 0.5)
				}
			}
		}
	}

	Data$Event <- 1

	Data$X4 <- sample(1:5,n,replace = TRUE)
	Data$X5 <- sample(c(1,2),n,replace = TRUE)
	Data$X6 <- runif(n,0,2)
	return(Data)
}

