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
	Train$ID <- predict(rtree, type="vector")
	Keys <- unique(Train$ID)
	Keys.MM <- matrix(c(Keys,1:length(Keys)),ncol=2)

	List.KM <- list()
	List.Med <- list()

	for(p in Keys){
		subset <- Train[Train$ID == p,]
		KM <-  survfit(Formula, data = subset)
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
source("LTRC_structure_grt.r")
source("Hothorn_gnrt.r")


Pred_LTRC <- function(N=300, distribution = "Exponential", C.rate=1, Trunc.max = 2, model = 1){

	set.seed(3)

	BS.ctree <- rep(NA,100)
	BS.ctree.RC <- rep(NA,100)
	BS.rpart <- rep(NA,100)
	BS.rpart.RC <- rep(NA,100)
	BScox <- rep(NA,100)
	BScox.RC <- rep(NA,100)

	Med.ctree <- rep(NA,100)
	Med.ctree.RC <- rep(NA,100)
	Med.rpart <- rep(NA,100)
	Med.rpart.RC <- rep(NA,100)
	Medcox <- rep(NA,100)
	Medcox.RC <- rep(NA,100)

	U.ctree <- rep(NA,100)
	U.ctree.RC <- rep(NA,100)
	U.rpart <- rep(NA,100)
	U.rpart.RC <- rep(NA,100)
	Ucox <- rep(NA,100)
      Ucox.RC <- rep(NA,100)

	for(W in 1:100){

		if(model == 1){
			Data.train <- LTRC.generate(n = N, Dist = distribution, censor.type = C.rate, truncation = Trunc.max)
			Data.test <- LTRC.generate.test(n=N, Dist = distribution)
		}else{
			Data.train <- Hothorn_gnrt(N, Model.type = model, Dist = distribution, censor.type=C.rate, truncation = Trunc.max)
			Data.test <- Hothorn_gnrt_test(N, Model.type = model,  Dist = distribution)
		}

		Ctree <- LTRC_ctree(Surv(Start,Obs,Event)~X1+X2+X3+X4+X5+X6, Data.train)
	      Ctree.RC <- ctree(Surv(Obs,Event)~X1+X2+X3+X4+X5+X6, Data.train)
            Rtree <- rpart.LTRC(Surv(Start,Obs,Event)~X1+X2+X3+X4+X5+X6, Data.train)
	      Rtree.RC <- rpart(Surv(Obs,Event)~X1+X2+X3+X4+X5+X6, Data.train,control = rpart.control(cp = 0.001))
		##==============
		if( length(unique(Rtree.RC$where))!=1 ){
			cventry.RC <- which.min(Rtree.RC$cptable[, "xerror"])
			cpcv.RC <- Rtree.RC$cptable[cventry.RC, "CP"]
			Rtree.RC <- prune(Rtree.RC,cp=cpcv.RC)
		}
		##==============
		Coxfit <- coxph(Surv(Start, Obs, Event)~X1+X2+X3+X4+X5+X6,data = Data.train)
		Coxfit.RC <- coxph(Surv(Obs, Event)~X1+X2+X3+X4+X5+X6,data = Data.train)

		km.LTRC0 <- survfit(Surv(Start, Obs, Event)~1, data = Data.train)
		km.RC0 <- survfit(Surv(Obs, Event)~1, data = Data.train)
			
		
		########################---L1 loss----#########################
		###------------------LTRC Ctree
		if(length(Ctree)==1){
			Med.ctree[W]<- mean(abs(rep(read.table(textConnection(capture.output(km.LTRC0)),skip=2,header=TRUE)$median,nrow(Data.test))-Data.test$Obs))
		}else{
			Med.ctree[W] <- mean(abs(unname(predict(Ctree, newdata=Data.test,type = "response"))-Data.test$Obs))
		}
		Null.LTRC <- mean(abs(rep(read.table(textConnection(capture.output(km.LTRC0)),skip=2,header=TRUE)$median,nrow(Data.test))-Data.test$Obs))
		U.ctree[W] <- 1 - Med.ctree[W]/Null.LTRC

		###------------------Ctree.RC
		if(length(Ctree.RC)==1){
			Med.ctree.RC[W]<- mean(abs(rep(read.table(textConnection(capture.output(km.RC0)),skip=2,header=TRUE)$median,nrow(Data.test))-Data.test$Obs))
		}else{
			Med.ctree.RC[W] <- mean(abs(unname(predict(Ctree.RC, newdata=Data.test,type = "response"))-Data.test$Obs))
		}
		Null.RC <- mean(abs(rep(read.table(textConnection(capture.output(km.RC0)),skip=2,header=TRUE)$median,nrow(Data.test))-Data.test$Obs))
		U.ctree.RC[W] <- 1-Med.ctree.RC[W]/Null.RC

		###------------------LTRC Rtree
		if(length(unique(Rtree$where))==1){
			Med.rpart[W] <- mean(abs(rep(read.table(textConnection(capture.output(km.LTRC0)),skip=2,header=TRUE)$median,nrow(Data.test))-Data.test$Obs))
		}else{
			Test.result <- Pred.rpart(Surv(Start,Obs,Event)~X1+X2+X3+X4+X5+X6, Data.train, Data.test)
			Med.rpart[W] <- mean(abs(Test.result$Medians-Data.test$Obs))
		}
		U.rpart[W] <- 1-Med.rpart[W]/Null.LTRC

		###------------------LTRC Rtree.RC
		if(length(unique(Rtree.RC$where))==1){
			Med.rpart.RC[W] <- mean(abs(rep(read.table(textConnection(capture.output(km.RC0)),skip=2,header=TRUE)$median,nrow(Data.test))-Data.test$Obs))
		}else{
			Test.result.RC <- Pred.rpart(Surv(Obs,Event)~X1+X2+X3+X4+X5+X6, Data.train, Data.test)
			Med.rpart.RC[W] <- mean(abs(Test.result.RC$Medians-Data.test$Obs))
		}
		U.rpart.RC[W] <- 1 - Med.rpart.RC[W]/Null.RC

		##-------------------LTRC Cox
		Medcox[W] <- mean(abs(unname(summary(survfit(Coxfit ,newdata=Data.test))$table[,5])-Data.test$Obs))
		Ucox[W] <- 1-Medcox[W]/Null.LTRC
	
		##-------------------Cox 
		Medcox.RC[W] <- mean(abs(unname(summary(survfit(Coxfit.RC ,newdata=Data.test))$table[,5])-Data.test$Obs))
		Ucox.RC[W] <- 1-Medcox.RC[W]/Null.RC

		#####################---Brier Score for Ctree----###########################
		if(length(Ctree)==1){
			Test.obj <- Surv(Data.test$Obs, Data.test$Event)
			BS.ctree[W] <- unname(sbrier(Test.obj, km.LTRC0)[1])
		}else{
			Test.obj <- Surv(Data.test$Obs, Data.test$Event)
			Pred.LTRC <- predict(Ctree, newdata=Data.test, type = "prob")
			BS.ctree[W]  <- unname(sbrier(Test.obj, Pred.LTRC)[1])
		}
		#####################---Brier Score for Ctree.RC----###########################
		if(length(Ctree.RC)==1){
			Test.obj <- Surv(Data.test$Obs, Data.test$Event)
			BS.ctree.RC[W] <- unname(sbrier(Test.obj, km.RC0)[1])
		}else{
			Test.obj <- Surv(Data.test$Obs, Data.test$Event)
			Pred.RC <- predict(Ctree.RC, newdata=Data.test, type = "prob")
			BS.ctree.RC[W]  <- unname(sbrier(Test.obj, Pred.RC)[1])
		}
		#####################---Brier Score for Rtree----###########################
		Test.obj <- Surv(Data.test$Obs, Data.test$Event)

		if(length(unique(Rtree$where))==1){
			BS.rpart[W] <- unname(sbrier(Test.obj, km.LTRC0)[1])
		}else{
			Test.KM <- Pred.rpart(Surv(Start,Obs,Event)~X1+X2+X3+X4+X5+X6, Data.train, Data.test)
			BS.rpart[W]  <- unname(sbrier(Test.obj, Test.KM$KMcurvs)[1])
		}
		#####################---Brier Score for Rtree.RC----###########################
		if(length(unique(Rtree.RC$where))==1){
			BS.rpart.RC[W] <- unname(sbrier(Test.obj, km.RC0)[1])
		}else{
			Test.KM.RC <- Pred.rpart(Surv(Obs,Event)~X1+X2+X3+X4+X5+X6, Data.train, Data.test)
			BS.rpart.RC[W]  <- unname(sbrier(Test.obj, Test.KM.RC$KMcurvs)[1])
		}
		########################---Brier Score for Coxph----#########################
		KM1 <- survfit(Coxfit, newdata = Data.test)

		Strat1 <- list()
		for( i in 1:nrow(Data.test) ){
			Strat1 <- c(Strat1, list(KM1[i]) )
		}

		Test.obj <- Surv(Data.test$Obs, Data.test$Event)
		BScox[W] <- unname(sbrier(Test.obj, Strat1)[1])
		########################---Brier Score for Coxph.RC----#########################
		KM1.RC <- survfit(Coxfit.RC, newdata = Data.test)

		Strat1.RC <- list()
		for( i in 1:nrow(Data.test) ){
			Strat1.RC <- c(Strat1.RC, list(KM1.RC[i]) )
		}

		Test.obj <- Surv(Data.test$Obs, Data.test$Event)
		BScox.RC[W] <- unname(sbrier(Test.obj, Strat1.RC)[1])
		#####################---Brier Score for Tree.RC----###########################
	}##End of for loop

	Brier_score <- cbind(BS.ctree, BS.ctree.RC, BS.rpart, BS.rpart.RC, BScox, BScox.RC)
	L1_loss <-  cbind(Med.ctree, Med.ctree.RC, Med.rpart, Med.rpart.RC, Medcox, Medcox.RC)
  	U_score <- cbind(U.ctree, U.ctree.RC, U.rpart, U.rpart.RC, Ucox, Ucox.RC)

	if( length(which(L1_loss == Inf))+length(which(L1_loss == -Inf)) != 0 ){
		ID1 <- which(L1_loss == Inf)
		ID2 <- which(L1_loss == -Inf)
		ID <- c(ID1,ID2)
		L1_loss[ID] <- NA
	}

	if( length(which(U_score == Inf))+length(which(U_score == -Inf)) != 0 ){
		ID3 <- which(U_score == Inf)
		ID4 <- which(U_score == -Inf)
		ID <- c(ID3,ID4)
		U_score[ID] <- NA
	}

  	BS.p1 <- wilcox.test( Brier_score[,1],  Brier_score[,3], paired = TRUE, alternative = "less")$p.value  ##ctree vs rpart
  	BS.p2 <- wilcox.test( Brier_score[,1],  Brier_score[,5], paired = TRUE, alternative = "less")$p.value  ##ctree vs cox
  	BS.p3 <- wilcox.test( Brier_score[,3],  Brier_score[,5], paired = TRUE, alternative = "less")$p.value  ##rpart vs cox
	BS.p4 <- wilcox.test( Brier_score[,1],  Brier_score[,2], paired = TRUE, alternative = "less")$p.value  ##ctree two version
  	BS.p5 <- wilcox.test( Brier_score[,3],  Brier_score[,4], paired = TRUE, alternative = "less")$p.value  ##rpart two version
  	BS.p6 <- wilcox.test( Brier_score[,5],  Brier_score[,6], paired = TRUE, alternative = "less")$p.value  ##cox two version
	BS.pv <- c(BS.p1,BS.p2,BS.p3,BS.p4,BS.p5,BS.p6)

	Med.p1 <- wilcox.test( L1_loss[,1],  L1_loss[,3], paired = TRUE, alternative = "less")$p.value  ##ctree vs rpart
  	Med.p2 <- wilcox.test( L1_loss[,1],  L1_loss[,5], paired = TRUE, alternative = "less")$p.value  ##ctree vs cox
  	Med.p3 <- wilcox.test( L1_loss[,3],  L1_loss[,5], paired = TRUE, alternative = "less")$p.value  ##rpart vs cox
	Med.p4 <- wilcox.test( L1_loss[,1],  L1_loss[,2], paired = TRUE, alternative = "less")$p.value  ##ctree two version
  	Med.p5 <- wilcox.test( L1_loss[,3],  L1_loss[,4], paired = TRUE, alternative = "less")$p.value  ##rpart two version
  	Med.p6 <- wilcox.test( L1_loss[,5],  L1_loss[,6], paired = TRUE, alternative = "less")$p.value  ##cox two version
	Med.pv <- c(Med.p1,Med.p2,Med.p3,Med.p4,Med.p5,Med.p6)

	U.p1 <- wilcox.test( U_score[,1],  U_score[,3], paired = TRUE, alternative = "less")$p.value  ##ctree vs rpart
  	U.p2 <- wilcox.test( U_score[,1],  U_score[,5], paired = TRUE, alternative = "less")$p.value  ##ctree vs cox
  	U.p3 <- wilcox.test( U_score[,3],  U_score[,5], paired = TRUE, alternative = "less")$p.value  ##rpart vs cox
	U.p4 <- wilcox.test( U_score[,1],  U_score[,2], paired = TRUE, alternative = "less")$p.value  ##ctree two version
  	U.p5 <- wilcox.test( U_score[,3],  U_score[,4], paired = TRUE, alternative = "less")$p.value  ##rpart two version
  	U.p6 <- wilcox.test( U_score[,5],  U_score[,6], paired = TRUE, alternative = "less")$p.value  ##cox two version
	U.pv <- c(U.p1,U.p2,U.p3,U.p4,U.p5,U.p6)
	
  	return(list(BS = Brier_score, L1= L1_loss, U = U_score, BS.PV = BS.pv, Med.PV = Med.pv, U.PV = U.pv))
}

###================================================================================================
###=================================================================================================
###model = 1 means data comes from tree structure, i.e. setting 1
###model = 2 means data comes from PH model, i.e. setting 2
###model = 3 means data comes from complex model, i.e. setting 3
###=================================================================================================

WebD.1.1 <- Pred_LTRC(N = 300, distribution = "Weibull-D", C.rate = 1, Trunc.max= 1, model = 1)
WebD.2.1 <- Pred_LTRC(N = 300, distribution = "Weibull-D", C.rate = 2, Trunc.max = 1, model = 1)
WebD.1.2 <- Pred_LTRC(N = 300, distribution = "Weibull-D", C.rate= 1, Trunc.max= 2, model = 1)
WebD.2.2 <- Pred_LTRC(N = 300, distribution = "Weibull-D", C.rate = 2, Trunc.max = 2, model = 1)
WebD.1.3 <- Pred_LTRC(N = 300, distribution = "Weibull-D", C.rate = 1, Trunc.max = 3, model = 1)
WebD.2.3 <- Pred_LTRC(N = 300, distribution = "Weibull-D", C.rate= 2, Trunc.max = 3, model = 1)

Bath.1.1 <- Pred_LTRC(N = 300, distribution = "Bathtub", C.rate = 1, Trunc.max = 1, model = 1)
Bath.2.1 <- Pred_LTRC(N = 300, distribution = "Bathtub", C.rate = 2, Trunc.max = 1, model = 1)
Bath.1.2 <- Pred_LTRC(N = 300, distribution = "Bathtub", C.rate = 1, Trunc.max = 2, model = 1)
Bath.2.2 <- Pred_LTRC(N = 300, distribution = "Bathtub", C.rate = 2, Trunc.max= 2, model = 1)
Bath.1.3 <- Pred_LTRC(N = 300, distribution = "Bathtub", C.rate = 1, Trunc.max= 3, model = 1)
Bath.2.3 <- Pred_LTRC(N = 300, distribution = "Bathtub", C.rate = 2, Trunc.max= 3, model = 1)

Exp.1.1 <- Pred_LTRC(N = 300, distribution = "Exponential", C.rate = 1, Trunc.max= 1, model = 1)
Exp.2.1 <- Pred_LTRC(N = 300, distribution = "Exponential", C.rate = 2, Trunc.max = 1, model = 1)
Exp.1.2 <- Pred_LTRC(N = 300, distribution = "Exponential", C.rate = 1, Trunc.max = 2, model = 1)
Exp.2.2 <- Pred_LTRC(N = 300, distribution = "Exponential", C.rate = 2, Trunc.max= 2, model = 1)
Exp.1.3 <- Pred_LTRC(N = 300, distribution = "Exponential", C.rate = 1, Trunc.max= 3, model = 1)
Exp.2.3 <- Pred_LTRC(N = 300, distribution = "Exponential", C.rate= 2, Trunc.max = 3, model = 1)

WebI.1.1 <- Pred_LTRC(N = 300, distribution = "Weibull-I", C.rate = 1, Trunc.max= 1, model = 1)
WebI.2.1 <- Pred_LTRC(N = 300, distribution = "Weibull-I", C.rate = 2, Trunc.max= 1, model = 1)
WebI.1.2 <- Pred_LTRC(N = 300, distribution = "Weibull-I", C.rate = 1, Trunc.max= 2, model = 1)
WebI.2.2 <- Pred_LTRC(N = 300, distribution = "Weibull-I", C.rate = 2, Trunc.max= 2, model = 1)
WebI.1.3 <- Pred_LTRC(N = 300, distribution = "Weibull-I", C.rate = 1, Trunc.max= 3, model = 1)
WebI.2.3 <- Pred_LTRC(N = 300, distribution = "Weibull-I", C.rate = 2, Trunc.max = 3, model = 1)

Log.1.1 <- Pred_LTRC(N = 300, distribution = "Lognormal", C.rate = 1, Trunc.max= 1, model = 1)
Log.2.1 <- Pred_LTRC(N = 300, distribution = "Lognormal", C.rate= 2, Trunc.max= 1, model = 1)
Log.1.2 <- Pred_LTRC(N = 300, distribution = "Lognormal", C.rate= 1, Trunc.max= 2, model = 1)
Log.2.2 <- Pred_LTRC(N = 300, distribution = "Lognormal", C.rate= 2, Trunc.max = 2, model = 1)
Log.1.3 <- Pred_LTRC(N = 300, distribution = "Lognormal", C.rate = 1, Trunc.max = 3, model = 1)
Log.2.3 <- Pred_LTRC(N = 300, distribution = "Lognormal", C.rate= 2, Trunc.max = 3, model = 1)


###==========================PLOT Fig.3 in paper=======================================

pdf("BS_light_300.pdf", width=10)


par(mfrow=c(3,5))
boxplot.matrix(Exp.1.1$BS, names=c("1","2","3","4","5","6"),main="Exponential")
boxplot.matrix(WebI.1.1$BS,names=c("1","2","3","4","5","6"), main="Weibull-I")
boxplot.matrix(WebD.1.1$BS, names=c("1","2","3","4","5","6"),main="Weibull-D")
boxplot.matrix(Log.1.1$BS, names=c("1","2","3","4","5","6"),main="Lognormal")
boxplot.matrix(Bath.1.1$BS, names=c("1","2","3","4","5","6"),main="Bathtub")
boxplot.matrix(Exp.1.2$BS, names=c("1","2","3","4","5","6"))
boxplot.matrix(WebI.1.2$BS, names=c("1","2","3","4","5","6"))
boxplot.matrix(WebD.1.2$BS, names=c("1","2","3","4","5","6"))
boxplot.matrix(Log.1.2$BS, names=c("1","2","3","4","5","6"))
boxplot.matrix(Bath.1.2$BS, names=c("1","2","3","4","5","6"))
boxplot.matrix(Exp.1.3$BS, names=c("1","2","3","4","5","6"))
boxplot.matrix(WebI.1.3$BS, names=c("1","2","3","4","5","6"))
boxplot.matrix(WebD.1.3$BS, names=c("1","2","3","4","5","6"))
boxplot.matrix(Log.1.3$BS, names=c("1","2","3","4","5","6"))
boxplot.matrix(Bath.1.3$BS, names=c("1","2","3","4","5","6"))


dev.off()




