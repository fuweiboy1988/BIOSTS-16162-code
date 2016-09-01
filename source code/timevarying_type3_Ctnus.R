
## Range_T function returns simulated survival time 
Range_T_ctn <- function(T1, T2, T3, DIST, X, U){
	x1 = X
	u = U
	if( !(T1 < T2 && T2 < T3) ){
			TTT <- sort(c(T1,T2,T3)) 
			t1 = TTT[1]
			t2 = TTT[2]
			t3 = TTT[3]
	}else{
		t1 = T1;
		t2 = T2;
		t3 = T3
	}

	if(DIST == "Exponential"){
		Lamda = 0.1
		Beta1 <- 0.8
  		Beta2 <- 1.4
		
		R1 = Lamda*exp(Beta1*x1)*t1
		R2 = Lamda*exp(Beta1*x1)*(t1+exp(Beta2)*(t2-t1))
		R3 = Lamda*exp(Beta1*x1)*(t1+exp(Beta2)*(t2-t1)+t3-t2)
		VEC = c(0,R1,R2,R3,Inf)

		R.ID <- findInterval(-log(u), VEC)
		if(R.ID == 1){
			Time = -log(u)/(Lamda*exp(Beta1*x1))
		}else if(R.ID == 2){
			Time = (-log(u)-Lamda*exp(Beta1*x1)*t1+Lamda*exp(Beta1*x1+Beta2)*t1)/(Lamda*exp(Beta1*x1+Beta2))
		}else if(R.ID == 3){
			Time = (-log(u)-Lamda*exp(Beta1*x1)*t1-Lamda*exp(Beta1*x1+Beta2)*(t2-t1)+Lamda*exp(Beta1*x1)*t2)/(Lamda*exp(Beta1*x1))
		}else if(R.ID == 4){
			Time = (-log(u)-Lamda*exp(Beta1*x1)*t1-Lamda*exp(Beta1*x1+Beta2)*(t2-t1)-Lamda*exp(Beta1*x1)*(t3-t2)+Lamda*exp(Beta1*x1+Beta2)*t3)/(Lamda*exp(Beta1*x1+Beta2))
		}else{
			stop("Cannot find interval where u lays in")
		}

	}else if(DIST == "Weibull"){
		Lamda = 0.3
		Beta1 = 0.9
  		Beta2 = 1.6
		V = 0.8

		R1 = Lamda*exp(Beta1*x1)*t1^V
		R2 = Lamda*exp(Beta1*x1)*(t1^V+exp(Beta2)*(t2^V-t1^V))
		R3 = Lamda*exp(Beta1*x1)*(t1^V+exp(Beta2)*(t2^V-t1^V)+t3^V-t2^V)
		VEC = c(0,R1,R2,R3,Inf)

		R.ID <- findInterval(-log(u), VEC)
		if(R.ID == 1){
			TT = -log(u)/(Lamda*exp(Beta1*x1))
		}else if(R.ID == 2){
			TT = (-log(u)-Lamda*exp(Beta1*x1)*t1^V+Lamda*exp(Beta1*x1+Beta2)*t1^V)/(Lamda*exp(Beta1*x1+Beta2))
		}else if(R.ID == 3){
			TT = (-log(u)-Lamda*exp(Beta1*x1)*t1^V-Lamda*exp(Beta1*x1+Beta2)*(t2^V-t1^V)+Lamda*exp(Beta1*x1)*t2^V)/(Lamda*exp(Beta1*x1))
		}else if(R.ID == 4){
			TT = (-log(u)-Lamda*exp(Beta1*x1)*t1^V-Lamda*exp(Beta1*x1+Beta2)*(t2^V-t1^V)-Lamda*exp(Beta1*x1)*(t3^V-t2^V)+Lamda*exp(Beta1*x1+Beta2)*t3^V)/(Lamda*exp(Beta1*x1+Beta2))
		}else{
			stop("Cannot find interval where u lays in")
		}

		Time = TT^(1/V)

	}else if(DIST == "Gompertz"){
		Lamda = 0.2
		Alpha = 0.1
		Beta1 = 1.2
  		Beta2 = 2.0
		
		R1 = Lamda*exp(Beta1*x1)/Alpha*(exp(Alpha*t1)-1)
		R2 = Lamda*exp(Beta1*x1)/Alpha*(exp(Alpha*t1)-1+exp(Beta2+Alpha*t2)-exp(Beta2+Alpha*t1))
		R3 = Lamda*exp(Beta1*x1)/Alpha*(exp(Alpha*t1)-1+exp(Beta2+Alpha*t2)-exp(Beta2+Alpha*t1)+exp(Alpha*t3)-exp(Alpha*t2))
		VEC = c(0,R1,R2,R3,Inf)

		R.ID <- findInterval(-log(u), VEC)
		if(R.ID == 1){
			T_T = 1+Alpha*(-log(u))/(Lamda*exp(Beta1*x1))
		}else if(R.ID == 2){
			T_T = Alpha*(-log(u))/(Lamda*exp(Beta1*x1+Beta2)) - (exp(Alpha*t1)-1-exp(Beta2+Alpha*t1))/exp(Beta2)
		}else if(R.ID == 3){
			T_T = 1+Alpha*(-log(u))/(Lamda*exp(Beta1*x1))-exp(Alpha*t1)-exp(Beta2+Alpha*t2)+exp(Beta2+Alpha*t1)+exp(Alpha*t2)
		}else if(R.ID == 4){
			T_T = Alpha*(-log(u))/(Lamda*exp(Beta1*x1+Beta2)) - (exp(Alpha*t1)-1+exp(Beta2+Alpha*t2)-exp(Beta2+Alpha*t1)+exp(Alpha*t3)-exp(Alpha*t2)-exp(Beta2+Alpha*t3))/exp(Beta2)
		}else{
			stop("Cannot find interval where u lays in")
		}

		Time = 1/Alpha*log(T_T)
	}

	result = list(T = Time, Row = R.ID)
	return(result)
}

###================================================================
Timevarying_gnrt3_ctn <- function(N = 200, Distribution = "Exponential", censor.rate = 1){
	Data <- as.data.frame(matrix(NA,4*N,9))
	names(Data)<-c("ID","X1","X2","X3","X4","X5","Start","Stop","Event")
	Data$ID <- rep(1:N,each=4)
	Data$X1 <- rep(sample(c(0,1),N,replace=TRUE),each=4)
	Data$X2 <- as.vector(rbind(runif(2*N,0,5), runif(2*N,5,10)))
	Data$X3 <- rep(sample(c(0,1),N,replace=TRUE),each=4)
	Data$X4 <- runif(4*N)	
	Data$X5 <- sample(1:5,4*N,replace=TRUE)

	Count = 1
	
	while(Count <= N){

		Ts = runif(3,0.6,6)

		TS = sort(Ts)
		t1 = TS[1]
		t2 = TS[2]
		t3 = TS[3]
		u = runif(1)
		x = unique(Data[Data$ID==Count,]$X1) 
		Data[Data$ID==Count,]$Start <- c(0,t1,t2,t3)
		Data[Data$ID==Count,]$Stop <- c(t1,t2,t3,NA)

		RT <- Range_T_ctn(t1,t2,t3,Distribution,x,u)
		Time = RT$T
		Row.ID = RT$Row
		
		if(Row.ID==1){
			Data[Data$ID==Count,][1,]$Stop = Time
			Data[Data$ID==Count,][1,]$Event = 1
		}else{
			Data[Data$ID==Count,][1:(Row.ID-1),]$Event=0
			Data[Data$ID==Count,][Row.ID,]$Event = 1
			Data[Data$ID==Count,][Row.ID,]$Stop = Time
		}
		
		Count = Count+1
	}##end of while loop

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
			Censor.time = rexp(N,rate = 1/15)
		}else if(Distribution == "Weibull"){
			Censor.time = rexp(N,rate = 1/7.8)
		}else if(Distribution == "Gompertz"){
			Censor.time = rexp(N,rate = 1/7.1)
		}
	}else if(censor.rate == 2){
		if(Distribution == "Exponential"){
			Censor.time = rexp(N,rate = 1/4.2)
		}else if(Distribution == "Weibull"){
			Censor.time = rexp(N,rate = 1/2.0)
		}else if(Distribution == "Gompertz"){
			Censor.time = rexp(N,rate = 1/1.8)
		}
	}else{
		stop("Wrong censoring type")
	}

	for( j in 1:length(unique(DATA$ID)) ){
		Vec <- c(0,DATA[DATA$ID==j,]$Stop,Inf)
		ID <- findInterval(Censor.time[j], Vec)
		
		if( ID <= nrow(DATA[DATA$ID==j,]) ){
			DATA[DATA$ID==j,][ID,]$C = 1
			DATA[DATA$ID==j,][ID,]$Event = 0
			DATA[DATA$ID==j,][ID,]$Stop = Censor.time[j]
			if( ID != nrow(DATA[DATA$ID==j,]) ){
				DATA[DATA$ID==j,][(ID+1):nrow(DATA[DATA$ID==j,]),]$Event = NA
			}
		}
	}

	Data1 <- DATA[!is.na(DATA$Event),]

	if(length(unique(Data1$ID))!=N){
		stop("ID length NOT equal to N")
	}

	return(Data1)
}##end of function

##======================================================================================

Timevarying_gnrt3_ctn_test <- function(N = 200, Distribution = "Exponential"){
	Data <- as.data.frame(matrix(NA,N,7))
	names(Data)<-c("X1","X2","X3","X4","X5","Obs","Event")
	Data$X1 <- sample(c(0,1),N,replace=TRUE)
	Data$X2 <- runif(N,0,10)
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
			
			if( Data[count,"X2"]>5 ){
				Coef = Data[count,"X1"]*Beta1 + (Data[count,"X2"]>5)*Beta2
				t = -log(u)/(Lamda*exp(Coef))
				if(t>0.6){ ##left truncated at t=0.6
					Data[count,"Obs"] = t
					count = count+1
				}
			}else{
				Coef = Data[count,"X1"]*Beta1 + (Data[count,"X2"]>5)*Beta2
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
			if( Data[count,"X2"]>5 ){
				Coef = Data[count,"X1"]*Beta1 + (Data[count,"X2"]>5)*Beta2
				t = (-log(u)/(Lamda*exp(Coef)))^(1/V)
				if(t>0.6){ ##left truncated at t=0.6
					Data[count,"Obs"] = t
					count = count+1
				}
			}else{
				Coef = Data[count,"X1"]*Beta1 + (Data[count,"X2"]>5)*Beta2
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
			if( Data[count,"X2"]>5 ){
				Coef = Data[count,"X1"]*Beta1 + (Data[count,"X2"]>5)*Beta2
				t = 1/Alpha*log(1 - Alpha*log(u)/(Lamda*exp(Coef)))
				if(t>0.6){ ##left truncated at t=0.6
					Data[count,"Obs"] = t
					count = count+1
				}
			}else{
				Coef = Data[count,"X1"]*Beta1 + (Data[count,"X2"]>5)*Beta2
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
