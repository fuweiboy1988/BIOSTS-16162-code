bias.test.gnrt <- function(n=200, Dist = "Exponential", Ctype = 1){

	DATA1 <- as.data.frame(matrix(NA,n,5))
	names(DATA1)<-c("Start","T","C","Obs","Event")
	Count = 0
	Seq <- seq(0,1,0.1)

	DATA1$X1 <- runif(n,0,1)
	DATA1$X2 <- runif(n,0,1) ## original runif(n,0,2)
	DATA1$X3 <- sample(Seq,n,replace=TRUE)
	DATA1$X4 <- sample(c(0,1),n,replace=TRUE) ## originally sample(c(0,0.58),n,replace=TRUE)
	DATA1$X5 <- sample(c(0,1),n,replace=TRUE)

	while(Count < n){
		L <- runif(1,0,2)
		
		if(Dist == "Exponential"){
			t <- rexp(1,0.25)
		}else if(Dist == "Weibull"){
			t <- rweibull(1, 0.9, 2)
		}else if(Dist == "Lognormal"){
			t <- rlnorm(1, meanlog = 1.4, sdlog = 0.4)
		}else{ stop("Wrong distribution type")}
		
		if(t > L){
			Count = Count + 1
			DATA1[Count,"Start"] <- L
			DATA1[Count,"T"] <- t
			
			if(Dist == "Exponential"){
				if(Ctype==1){
					DATA1[Count,"C"] <- L + rexp(1,1/16)
				}else if(Ctype==2){
					DATA1[Count,"C"] <- L + rexp(1,1/4)
				}
			}else if(Dist == "Weibull"){
				 if(Ctype==1){
					DATA1[Count,"C"] <- L + rexp(1,1/9)
				}else if(Ctype==2){
					DATA1[Count,"C"] <- L + rexp(1,1/2.4)
				}

			}else if(Dist == "Lognormal"){
				if(Ctype==1){
					DATA1[Count,"C"] <- L + rexp(1,1/14.3)
				}else if(Ctype==2){
					DATA1[Count,"C"] <- L + rexp(1,1/4.3)
				}
			}

			if( t <= DATA1[Count,"C"]){
				DATA1[Count,"Obs"] <- t
				DATA1[Count,"Event"] <- 1
			}else{
				DATA1[Count,"Obs"] <- DATA1[Count,"C"]
				DATA1[Count,"Event"] <- 0
			}
		}
	}#end of while	
	return(DATA1)
}