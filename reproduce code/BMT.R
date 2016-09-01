## make sure to set the directory correctly 
setwd("C:/Users/wfu/Desktop/biostatistics/code")

library(survival)
library(KMsurv)
source("LTRC_ctree.R")
source("rpart_LTRC.R")


##Loading the BMT data
data(bmt)

BMT <- as.data.frame(matrix(NA,1,18))

for(i in 1:(nrow(bmt))){
	X <- bmt[i,];
	X2 <- sort(unique(c(min(X$t2,X$ta),min(X$t2,X$tc),min(X$t2,X$tp))))
	if(length(X2)==1 && X2 == X$t2){
		BMT <- rbind(BMT,c(i,0,X$t2,X$d3,0,0,0,as.numeric(X[,c(1,13:22)])))
	}else{
		Cut.p <- unique(c(0,X2,X$t2))
		Obj <- matrix(NA,length(Cut.p)-1,18)
		for(k in 1:(length(Cut.p)-1)){
			Obj[k,8:18] = as.numeric(X[,c(1,13:22)])
		}
		Obj[,1] <- i;
		Obj[,2] <- Cut.p[1:(length(Cut.p)-1)];
		Obj[,3] <- Cut.p[2:length(Cut.p)];	
		Obj[,4] <- 0
		Obj[length(Cut.p)-1,4] = X$d3
		row.ta <- findInterval(X$ta,Cut.p)
		row.tc <- findInterval(X$tc,Cut.p)
		row.tp <- findInterval(X$tp,Cut.p)
		Obj[,5] <- 0
		if(row.ta <= (length(Cut.p)-1)){Obj[row.ta:(length(Cut.p)-1),5] =1}
		Obj[,6] <- 0
		if(row.tc <= (length(Cut.p)-1)){Obj[row.tc:(length(Cut.p)-1),6] =1}
		Obj[,7] <- 0
		if(row.tp <= (length(Cut.p)-1)){Obj[row.tp:(length(Cut.p)-1),7] =1}
		BMT<-rbind(BMT,Obj)
	}
}
BMT <- BMT[-1,]
names(BMT)<-c("ID","start","stop","event","aGVHD","cGVHD","Platelets", "Group", "z1","z2","z3","z4","z5","z6","z7","FAB","z9","z10")

BMT$FAB <- as.factor(BMT$FAB)
BMT$Group <- as.factor(BMT$Group)
BMT$Platelets <- as.factor(BMT$Platelets)
BMT$cGVHD <- as.factor(BMT$cGVHD)
BMT$aGVHD <- as.factor(BMT$aGVHD)


##Cox model result of BMT data
Coxfit.full <- coxph(Surv(start,stop, event) ~aGVHD+cGVHD+Platelets+FAB+Group, BMT[,-1])
Coxfit.full

###Cox model result on subset of BMT data, with dease group-AML Low Risk
coxph(Surv(start,stop, event) ~aGVHD+cGVHD+Platelets+FAB, BMT[BMT$Group==2,-1])

##test proportional assumption
cox.zph(Coxfit.full)

##plot the graphes of the Schoenfeld residuals against transformed time
par(mfrow=c(3, 2))
plot(cox.zph(Coxfit.full))


#####=========Tree results=======================##########

set.seed(1)
Ctree <-  LTRC_ctree(Surv(start,stop, event) ~ aGVHD+cGVHD+Platelets+FAB+Group, BMT[,-1])
plot(Ctree)

set.seed(12)
Rtree <- rpart.LTRC(Surv(start,stop, event)~ aGVHD+cGVHD+Platelets+FAB+Group, BMT[,-1])
###Convert the LTRCART object to compatible with party object
survobj <- Surv(BMT$start, BMT$stop, BMT$event)
py <- as.party(Rtree)
py$fitted[["(response)"]]<- survobj
plot(py)

