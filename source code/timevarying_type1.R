

Timevarying_gnrt <- function(N = 200, Distribution = "Exponential", censor.rate = 1){
	Data <- as.data.frame(matrix(NA,2*N,9))
	names(Data)<-c("ID","X1","X2","X3","X4","X5","Start","Stop","Event")
	Data$ID <- rep(1:N,each=2)
	Data$X1 <- rep(sample(c(0,1),N,replace=TRUE),each=2)
	Data$X2 <- rep(c(0,1),N)
	Data$X3 <- rep(sample(c(0,1),N,replace=TRUE),each=2)
	Data$X4 <- runif(2*N)	
	Data$X5 <- sample(1:5,2*N,replace=TRUE)

	Count = 1

	if(Distribution == "Exponential"){
		Lamda = 0.1
		Beta1 <- 0.8
  		Beta2 <- 1.4
		
		while(Count <= N){
			t0 = runif(1,0.6,6) ## sample t0
			u = runif(1,0,1) ## sample u
			x1 = unique(Data[Data$ID==Count,]$X1)

			if( -log(u) < Lamda*exp(Beta1*x1)*t0 ){
				Data[Data$ID==Count,][1,]$Start = 0
				T = -log(u)/(Lamda*exp(Beta1*x1))
				Data[Data$ID==Count,][1,]$Stop = T
				Data[Data$ID==Count,][1,]$Event = 1

			}else{
				Data[Data$ID==Count,][1,]$Start = 0
				Data[Data$ID==Count,][1,]$Stop = t0
				Data[Data$ID==Count,][1,]$Event = 0
				
				Data[Data$ID==Count,][2,]$Start = t0
				T = (-log(u)-Lamda*exp(Beta1*x1)*t0+Lamda*exp(Beta1*x1+Beta2)*t0)/(Lamda*exp(Beta1*x1+Beta2))
				Data[Data$ID==Count,][2,]$Stop = T
				Data[Data$ID==Count,][2,]$Event = 1
			}
			Count = Count+1
		}
		
	}else if(Distribution == "Weibull"){
		Lamda = 0.3
		Beta1 = 0.9
  		Beta2 = 1.6
		V = 0.8

		while(Count <= N){
			t0 = runif(1,0.6,6) ## runif(1,0.7,3.5)
			u = runif(1,0,1) ## sample u
			x1 = unique(Data[Data$ID==Count,]$X1)

			if( -log(u) < Lamda*exp(Beta1*x1)*t0^V ){
				Data[Data$ID==Count,][1,]$Start = 0
				T = (-log(u)/(Lamda*exp(Beta1*x1)))^(1/V)
				Data[Data$ID==Count,][1,]$Stop = T
				Data[Data$ID==Count,][1,]$Event = 1

			}else{
				Data[Data$ID==Count,][1,]$Start = 0
				Data[Data$ID==Count,][1,]$Stop = t0
				Data[Data$ID==Count,][1,]$Event = 0
				
				Data[Data$ID==Count,][2,]$Start = t0
				T1 = (-log(u)-Lamda*exp(Beta1*x1)*t0^V+Lamda*exp(Beta1*x1)*exp(Beta2)*t0^V)/(Lamda*exp(Beta1*x1)*exp(Beta2))
				T = T1^(1/V)
				Data[Data$ID==Count,][2,]$Stop = T
				Data[Data$ID==Count,][2,]$Event = 1
			}
			Count = Count+1
		}

	}else if(Distribution == "Gompertz"){
		Lamda = 0.2
		Alpha = 0.1
		Beta1 = 1.2
  		Beta2 = 2.0
		
		while(Count <= N){
			t0 = runif(1,0.6,6) ## sample t0
			u = runif(1,0,1) ## sample u
			x1 = unique(Data[Data$ID==Count,]$X1)

			if( -log(u) < Lamda*exp(Beta1*x1)/Alpha*(exp(Alpha*t0)-1) ){
				Data[Data$ID==Count,][1,]$Start = 0
				T1 = 1+Alpha*(-log(u))/(Lamda*exp(Beta1*x1))
				T = (1/Alpha)*log(T1)
				Data[Data$ID==Count,][1,]$Stop = T
				Data[Data$ID==Count,][1,]$Event = 1

			}else{
				Data[Data$ID==Count,][1,]$Start = 0
				Data[Data$ID==Count,][1,]$Stop = t0
				Data[Data$ID==Count,][1,]$Event = 0
				
				Data[Data$ID==Count,][2,]$Start = t0
				T1 = Alpha*(-log(u))/(Lamda*exp(Beta1*x1+Beta2))-(exp(Alpha*t0)-1-exp(Beta2+Alpha*t0))/exp(Beta2)
				T = (1/Alpha)*log(T1)
				Data[Data$ID==Count,][2,]$Stop = T
				Data[Data$ID==Count,][2,]$Event = 1
			}
			Count = Count+1
		}
	}

	DATA <- Data[!is.na(Data$Event),]
	###================== Add Censoring =========================================
	DATA$C <- 0

	if(length(unique(DATA$ID)) != N){
		stop("ID length NOT equal to N")
	}
	
	if(censor.rate == 0){
		Censor.time <- rep(Inf,N)
	}else if(censor.rate == 1){
		if(Distribution == "Exponential"){
			Censor.time = rexp(N,rate = 1/16)
		}else if(Distribution == "Weibull"){
			Censor.time = rexp(N,rate = 1/9)
		}else if(Distribution == "Gompertz"){
			Censor.time = rexp(N,rate = 1/8.6)
		}
	}else if(censor.rate == 2){
		if(Distribution == "Exponential"){
			Censor.time = rexp(N,rate = 1/4.5)
		}else if(Distribution == "Weibull"){
			Censor.time = rexp(N,rate = 1/2.2)
		}else if(Distribution == "Gompertz"){
			Censor.time = rexp(N,rate = 1/2.2)
		}
	}else{
		stop("Wrong censoring type")
	}
	
	for( j in 1:length(unique(DATA$ID)) ){
		if(nrow(DATA[DATA$ID==j,]) == 1){# one row case
			if( DATA[DATA$ID==j,]$Stop > Censor.time[j]){ ##Censored case
				DATA[DATA$ID==j,]$Stop = Censor.time[j]
				DATA[DATA$ID==j,]$Event = 0
				DATA[DATA$ID==j,]$C = 1
			}
		}else{##two rows case
			if( DATA[DATA$ID==j,][2,]$Stop > Censor.time[j] ){
				if(DATA[DATA$ID==j,][1,]$Stop > Censor.time[j]){
					DATA[DATA$ID==j,][2,]$Event = NA
					DATA[DATA$ID==j,][1,]$Stop = Censor.time[j]
					DATA[DATA$ID==j,][1,]$C = 1
				}else{
					DATA[DATA$ID==j,][2,]$Stop = Censor.time[j]
					DATA[DATA$ID==j,][2,]$Event = 0 
					DATA[DATA$ID==j,][2,]$C = 1 
				}
			}
		}
	}

	Data1 <- DATA[!is.na(DATA$Event),]
	if(length(unique(Data1$ID))!=N){
		stop("ID length NOT equal to N")
	}
	return(Data1)
}##end of function

###=======================================================================
Timevarying_gnrt_test <- function(N = 200, Distribution = "Exponential"){
	Data <- as.data.frame(matrix(NA,N,7))
	names(Data)<-c("X1","X2","X3","X4","X5","Obs","Event")
	Data$X1 <- sample(c(0,1),N,replace=TRUE)
	Data$X2 <- sample(c(0,1),N,replace=TRUE)
	Data$X3 <- sample(c(0,1),N,replace=TRUE)
	Data$X4 <- runif(N)	
	Data$X5 <- sample(1:5,N,replace=TRUE)
	Data$Event <- 1
	
	count=1

	while(count <= N){
		u = runif(1);

		if(Distribution == "Exponential"){
			Lamda = 0.1
			Beta1 <- 0.8
  			Beta2 <- 1.4
			
			if( Data[count,"X2"]>=1 ){
				Coef = Data[count,"X1"]*Beta1 + Data[count,"X2"]*Beta2
				t = -log(u)/(Lamda*exp(Coef))
				if(t>0.6){ ##left truncated at t=0.6
					Data[count,"Obs"] = t
					count = count+1
				}
			}else{
				Coef = Data[count,"X1"]*Beta1 + Data[count,"X2"]*Beta2
				t = -log(u)/(Lamda*exp(Coef))
				if(t>6){ ##right censored at t=6
					Data[count,"Obs"] = 6
					Data[count,"Event"] = 0
				}else{
					Data[count,"Obs"] = t
				}
				count=count+1
			}
			
		}else if(Distribution == "Weibull"){
			Lamda = 0.3
			Beta1 = 0.9
  			Beta2 = 1.6
			V = 0.8
			if( Data[count,"X2"]>=1 ){
				Coef = Data[count,"X1"]*Beta1 + Data[count,"X2"]*Beta2
				t = (-log(u)/(Lamda*exp(Coef)))^(1/V)
				if(t>0.6){ ##left truncated at t=0.6
					Data[count,"Obs"] = t
					count = count+1
				}
			}else{
				Coef = Data[count,"X1"]*Beta1 + Data[count,"X2"]*Beta2
				t = (-log(u)/(Lamda*exp(Coef)))^(1/V)
				if(t>6){ ##right censored at t=6
					Data[count,"Obs"] = 6
					Data[count,"Event"] = 0
				}else{
					Data[count,"Obs"] = t
				}
				count=count+1
			}
		}else if(Distribution == "Gompertz"){
			Lamda = 0.2
			Alpha = 0.1
			Beta1 = 1.2
  			Beta2 = 2.0
			if( Data[count,"X2"]>=1 ){
				Coef = Data[count,"X1"]*Beta1 + Data[count,"X2"]*Beta2
				t = 1/Alpha*log(1 - Alpha*log(u)/(Lamda*exp(Coef)))
				if(t>0.6){ ##left truncated at t=0.6
					Data[count,"Obs"] = t
					count = count+1
				}
			}else{
				Coef = Data[count,"X1"]*Beta1 + Data[count,"X2"]*Beta2
				t = 1/Alpha*log(1 - Alpha*log(u)/(Lamda*exp(Coef)))
				if(t>6){ ##right censored at t=6
					Data[count,"Obs"] = 6
					Data[count,"Event"] = 0
				}else{
					Data[count,"Obs"] = t
				}
				count=count+1
			}
		}else{stop("wrong distribution")}
	}

	return(Data)
}



