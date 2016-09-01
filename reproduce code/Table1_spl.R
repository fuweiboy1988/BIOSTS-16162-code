## make sure to set the directory correctly 
setwd("C:/Users/wfu/Desktop/biostatistics/code")

####=======================================================================
getUsefulPredictors <- function(x) {
  varid <- nodeapply(x, ids = nodeids(x), FUN = function(n) split_node(n)$varid)
  #varid <- unique(unlist(varid))
  varid <- unlist(varid)
  names(data_party(x))[varid]
}

getPositions <-function(x){
	varid <- nodeapply(x, ids = nodeids(x), FUN = function(n) split_node(n)$varid)
	varid <- unlist(varid)
      names(varid)
}
####=======================================================================
library(survival)
library(partykit)
library(rpart)
source("rpart_LTRC_noprune.R")
source("LTRC_ctree.R")
source("bias_test_grnt.R")

bias_test <- function(N=200, distribution = "Weibull", censor.rate = 1){

	set.seed(1)

	Ctree.1 = 0
	Ctree.2 = 0
	Ctree.3 = 0
	Ctree.4 = 0
	Ctree.5 = 0

	Rtree.1 = 0
	Rtree.2 = 0
	Rtree.3 = 0
	Rtree.4 = 0
	Rtree.5 = 0

	for(W in 1:10000){
		
		DATA1 <- bias.test.gnrt(n=N, Dist = distribution, Ctype = censor.rate)

	    	Ctree <- LTRC_ctree(Surv(Start,Obs, Event) ~ X1+X2+X3+X4+X5, DATA1, Control = ctree_control(mincriterion = 0,maxdepth=1))

		if( length(Ctree)!=3 || width(Ctree)!=2 ){
			stop("No root split returned for LTRC_ctree");
		}else{
			Varable = getUsefulPredictors(Ctree)

			if(Varable == "X1"){
				Ctree.1 = Ctree.1 + 1
			} else if(Varable == "X2"){
				Ctree.2 = Ctree.2 + 1
			} else if(Varable == "X3"){
				Ctree.3 = Ctree.3 + 1
			} else if(Varable == "X4"){
				Ctree.4 = Ctree.4 + 1
			} else if(Varable == "X5"){
				Ctree.5 = Ctree.5 + 1
			} 
		}

		Rtree <- rpart.LTRC(Surv(Start,Obs, Event) ~ X1+X2+X3+X4+X5, DATA1, control = rpart.control(cp=0,maxdepth=1))

		if( length(unique(Rtree$where))!=2 ){
			stop("No root split returned for LTRC_rpart");
		}else{
			if(Rtree$frame$var[1]=="X1"){
				Rtree.1 = Rtree.1 + 1;
			}else if(Rtree$frame$var[1]=="X2"){
				Rtree.2 = Rtree.2 + 1;
			}else if(Rtree$frame$var[1]=="X3"){
				Rtree.3 = Rtree.3 + 1;
			}else if(Rtree$frame$var[1]=="X4"){
				Rtree.4 = Rtree.4 + 1;
			}else if(Rtree$frame$var[1]=="X5"){
				Rtree.5 = Rtree.5 + 1;
			}
		}

	}##end of for loop
	
	result <- list( C.split = c(Ctree.1,Ctree.2,Ctree.3,Ctree.4,Ctree.5), R.split=c(Rtree.1,Rtree.2,Rtree.3,Rtree.4,Rtree.5))
	return(result)
}
#############################################################################
######## censor.rate = 1 means light censoring 
######## censor.rate = 2 means heavy censoring 


Exp1 <- bias_test(N = 200, distribution = "Exponential", censor.rate = 1)
Exp2 <- bias_test(N = 200, distribution = "Exponential", censor.rate = 2)
Log1 <- bias_test(N = 200, distribution = "Lognormal", censor.rate = 1)
Log2 <- bias_test(N = 200, distribution = "Lognormal", censor.rate = 2)
Web1 <- bias_test(N = 200, distribution = "Weibull", censor.rate = 1)
Web2 <- bias_test(N = 200, distribution = "Weibull", censor.rate = 2)

####====Comduct the Chi-sqaure test======
### C.split is the LTRCIT result
### R.split is the LTRCART result
####================================

chisq.test(Exp1$C.split,p = rep(1/5, 5))
chisq.test(Exp1$R.split,p = rep(1/5, 5))
chisq.test(Exp2$C.split,p = rep(1/5, 5))
chisq.test(Exp2$R.split,p = rep(1/5, 5))

chisq.test(Log1$C.split,p = rep(1/5, 5))
chisq.test(Log1$R.split,p = rep(1/5, 5))
chisq.test(Log2$C.split,p = rep(1/5, 5))
chisq.test(Log2$R.split,p = rep(1/5, 5))

chisq.test(Web1$C.split,p = rep(1/5, 5))
chisq.test(Web1$R.split,p = rep(1/5, 5))
chisq.test(Web2$C.split,p = rep(1/5, 5))
chisq.test(Web2$R.split,p = rep(1/5, 5))

