######################################################
## make sure to first load the necessary source code in the 
## "source code" file before running code below.
## If you are not sure which source code to load,
## load all of them just to be safe. 
######################################################

###================================================
Pred.rpart <- function(formula, train, test){

	if( length(formula[[2]])==3){
	      rtree <- rpart(formula, train)
		Formula = formula
		Formula[[3]] = 1

	}else if(length(formula[[2]])==4){
		rtree <- rpart.LTRC(formula, train)
		Formula = formula
		Formula[[3]] = 1

	}else{
		stop("Not Surv object")
	}
	
	Train = train
	Train$ID <- predict(rtree, type="vector") ##get the node ID for training data
	Keys <- unique(Train$ID)
	Keys.MM <- matrix(c(Keys,1:length(Keys)),ncol=2)

	List.KM <- list()
	List.Med <- list()

	for(p in Keys){
		subset <- Train[Train$ID == p,]
		KM <-  survfit(Formula, data = subset) ##Fit KM curve for each node in fitted tree
		Median <- read.table(textConnection(capture.output(KM)),skip=2,header=TRUE)$median
		List.KM[[ Keys.MM[Keys.MM[,1]==p,2] ]] = KM
		List.Med[[ Keys.MM[Keys.MM[,1]==p,2] ]] = Median
	}

	Test.keys <- predict(rtree, newdata = test, type="vector")
	Test.ID <- match(Test.keys,Keys.MM[,1])
	Test.KM <- List.KM[Test.ID]
	Test.Med <- unlist(List.Med[Test.ID])

	result <- list(KMcurvs = Test.KM, Medians = Test.Med)
	return(result)
}
###===============================================
library(pec)
library(ipred)
library(prodlim)
library(survAUC)
library(partykit)
library(rpart)
library(survival)

source("rpart_LTRC.R")
source("LTRC_ctree.R")
source("timevarying_type1.r")
source("timevarying_type3.r")
source("timevarying_type1_Ctnus.r")
source("timevarying_type3_Ctnus.r")

Pred_TV <- function(n=300, Dist = "Exponential", Jumptype = 1, Ctype = 1, Var = "Binary"){
	set.seed(1)
	
	BS.ctree <- rep(NA,500)
	BS.rpart <- rep(NA,500)	
	BScox <- rep(NA,500)

	W=1
	while(W <= 500){
		if(Var == "Binary"){
			if(Jumptype == 1){
				Data.train <- Timevarying_gnrt(N = n, Distribution = Dist, censor.rate = Ctype)
				Data.test <- Timevarying_gnrt_test(N=n, Distribution = Dist)
			}else if(Jumptype == 3){
				Data.train <- Timevarying_gnrt3(N = n, Distribution = Dist, censor.rate = Ctype)
				Data.test <- Timevarying_gnrt3_test(N=n, Distribution = Dist)
			}
		}else if(Var == "Continuous"){
			if(Jumptype == 1){
				Data.train <- Timevarying_gnrt_ctn(N = n, Distribution = Dist, censor.rate = Ctype)
				Data.test <- Timevarying_gnrt_ctn_test(N=n, Distribution = Dist)
			}else if(Jumptype == 3){
				Data.train <- Timevarying_gnrt3_ctn(N = n, Distribution = Dist, censor.rate = Ctype)
				Data.test <- Timevarying_gnrt3_ctn_test(N=n, Distribution = Dist)
			}
		}else{stop("Wrong split variable type")}

		
		Ctree <- LTRC_ctree(Surv(Start,Stop,Event)~X1+X2+X3+X4+X5, Data.train)
		Rtree <- rpart.LTRC(Surv(Start,Stop,Event)~X1+X2+X3+X4+X5, Data.train)
		Coxfit <- coxph(Surv(Start,Stop,Event)~X1+X2+X3+X4+X5, Data.train)
		
		km.LTRC0 <- survfit(Surv(Start,Stop,Event)~1, data = Data.train)
		#####################---Brier Score for Ctree----###########################
		if(length(Ctree)==1){
			Test.obj <- Surv(Data.test$Obs, Data.test$Event)
			BS.ctree[W] <- unname(sbrier(Test.obj, km.LTRC0)[1])
		}else{
			Test.obj <- Surv(Data.test$Obs, Data.test$Event)
			Pred.LTRC <- predict(Ctree, newdata=Data.test, type = "prob")
			BS.ctree[W]  <- unname(sbrier(Test.obj, Pred.LTRC)[1])
		}
		#####################---Brier Score for Rtree----###########################
		Test.obj <- Surv(Data.test$Obs, Data.test$Event)

		if(length(unique(Rtree$where))==1){
			BS.rpart[W] <- unname(sbrier(Test.obj, km.LTRC0)[1])
		}else{
			Test.KM <- Pred.rpart(Surv(Start,Stop,Event)~X1+X2+X3+X4+X5, Data.train, Data.test)
			BS.rpart[W]  <- unname(sbrier(Test.obj, Test.KM$KMcurvs)[1])
		}
		########################---Brier Score for Coxph----#########################
		KM1 <- survfit(Coxfit, newdata = Data.test)

		Strat1 <- list()
		for( i in 1:nrow(Data.test) ){
			Strat1 <- c(Strat1, list(KM1[i]) )
		}

		Test.obj <- Surv(Data.test$Obs, Data.test$Event)
		BScox[W] <- unname(sbrier(Test.obj, Strat1)[1])
		############################################################################
		W=W+1
	}##end of while loop
	Brier_score <- cbind(BS.ctree, BS.rpart, BScox)	
	
	BS.p1 <- wilcox.test( Brier_score[,1],  Brier_score[,2], paired = TRUE, alternative = "less")$p.value  ##ctree vs rpart
  	BS.p2 <- wilcox.test( Brier_score[,1],  Brier_score[,3], paired = TRUE, alternative = "less")$p.value  ##ctree vs cox
  	BS.p3 <- wilcox.test( Brier_score[,2],  Brier_score[,3], paired = TRUE, alternative = "less")$p.value  ##rpart vs cox
	BS.pv <- c(BS.p1,BS.p2,BS.p3)

	return(list(BS = Brier_score,BS.PV = BS.pv))
}	
###=======================================================================#####
###=======================================================================#####
### n is the sample size; Dist is the distribution; Jumptype = 1 for dichotomous
### type I and  Jumptype = 3 for type II; Ctype = 0 mean 0% censoring,
### 1 mean 20% censoring and 2 means 50% censoring; Var = "Binary"(default) means 
### time-varying variable X2 is binary, while "Continuous" means X2 is continuous
###=======================================================================#####
###=======================================================================#####
### Code below gives the IBS boxplot with binary X2, dichotomous I with n=100
### i.e.  Figure SM23
###=======================================================================#####
Exp.0.100 <- Pred_TV(n=100, Dist = "Exponential", Jumptype = 1, Ctype = 0)
Exp.1.100 <- Pred_TV(n=100, Dist = "Exponential", Jumptype = 1, Ctype = 1)
Exp.2.100 <- Pred_TV(n=100, Dist = "Exponential", Jumptype = 1, Ctype = 2)

Web.0.100 <- Pred_TV(n = 100, Dist = "Weibull", Jumptype = 1, Ctype = 0)
Web.1.100 <- Pred_TV(n = 100, Dist = "Weibull", Jumptype = 1, Ctype = 1)
Web.2.100 <- Pred_TV(n = 100, Dist = "Weibull", Jumptype = 1, Ctype = 2)

Gom.0.100 <- Pred_TV(n = 100, Dist = "Gompertz", Jumptype = 1, Ctype = 0)
Gom.1.100 <- Pred_TV(n = 100, Dist = "Gompertz", Jumptype = 1, Ctype = 1)
Gom.2.100 <- Pred_TV(n = 100, Dist = "Gompertz", Jumptype = 1, Ctype = 2)

###==========================PLOT=======================================

par(mfrow=c(3,3))
boxplot.matrix(Exp.0.100$BS, names=c("LTRCIT","LTRCART","Cox"), main="Exponential")
boxplot.matrix(Web.0.100$BS, names=c("LTRCIT","LTRCART","Cox"), main="Weibull")
boxplot.matrix(Gom.0.100$BS, names=c("LTRCIT","LTRCART","Cox"), main="Gompertz")
boxplot.matrix(Exp.1.100$BS, names=c("LTRCIT","LTRCART","Cox"))
boxplot.matrix(Web.1.100$BS, names=c("LTRCIT","LTRCART","Cox"))
boxplot.matrix(Gom.1.100$BS, names=c("LTRCIT","LTRCART","Cox"))
boxplot.matrix(Exp.2.100$BS, names=c("LTRCIT","LTRCART","Cox"))
boxplot.matrix(Web.2.100$BS, names=c("LTRCIT","LTRCART","Cox"))
boxplot.matrix(Gom.2.100$BS, names=c("LTRCIT","LTRCART","Cox"))


####=======Signed-rank test =========================

Exp.0.100$BS.PV
Web.0.100$BS.PV
Gom.0.100$BS.PV
Exp.1.100$BS.PV
Web.1.100$BS.PV
Gom.1.100$BS.PV
Exp.2.100$BS.PV
Web.2.100$BS.PV
Gom.2.100$BS.PV




