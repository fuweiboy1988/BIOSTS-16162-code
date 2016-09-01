## make sure to set the directory correctly 
setwd("C:/Users/wfu/Desktop/biostatistics/code")


library(survival)
library(KMsurv)
library(rpart.plot)
source("LTRC_ctree.R")
source("rpart_LTRC.R")

###=========================================================
### Baseline case, i.e. time-independent covariates
###=========================================================

pbc$status <- (pbc$status == 2)
pbc$albumin<-log(pbc$albumin)
pbc$bili<-log(pbc$bili)
pbc$protime<-log(pbc$protime)

coxph(Surv(time, status) ~ age+edema+alk.phos+chol+ast+platelet+spiders+hepato+ascites+albumin+bili+protime, pbc)

set.seed(1)
Ctree <- ctree(Surv(time, status) ~ age+edema+alk.phos+chol+ast+platelet+spiders+hepato+ascites+albumin+bili+protime, pbc)
plot(Ctree)

set.seed(13)
Rtree <- rpart(Surv(time, status) ~ age+edema+alk.phos+chol+ast+platelet+spiders+hepato+ascites+albumin+bili+protime, pbc)
## prune the rpart tree 
cventry <- which.min(Rtree$cptable[, "xerror"])
cpcv <- Rtree$cptable[cventry, "CP"]
Rtree <- prune(Rtree,cp=cpcv)
plot(as.party(Rtree))

###=========================================================
### Complete data case, i.e. time-varying covariates
###=========================================================
first <- with(pbcseq, c(TRUE, diff(id) !=0)) #first id for each subject
last  <- c(first[-1], TRUE)  #last id
time1 <- with(pbcseq, ifelse(first, 0, day))
time2 <- with(pbcseq, ifelse(last,  futime, c(day[-1], 0)))
event <- with(pbcseq, ifelse(last,  status, 0))
event <- 1*(event==2)

pbcseq$time1 <- time1
pbcseq$time2 <- time2
pbcseq$event <-  event

pbcseq$albumin<-log(pbcseq$albumin)
pbcseq$bili<-log(pbcseq$bili)
pbcseq$protime<-log(pbcseq$protime)

coxph(Surv(time1, time2, event) ~ age+edema+alk.phos+chol+ast+platelet+spiders+hepato+ascites+albumin+bili+protime, pbcseq)

set.seed(2)
Ctree <- LTRC_ctree(Surv(time1, time2, event) ~ age+edema+alk.phos+chol+ast+platelet+spiders+hepato+ascites+albumin+bili+protime, pbcseq)
plot(Ctree)

set.seed(2)
Rtree <- rpart.LTRC(Surv(time1, time2, event) ~ age+edema+alk.phos+chol+ast+platelet+spiders+hepato+ascites+albumin+bili+protime, pbcseq)
###Convert the LTRCART object to compatible with party object
survobj <- Surv(pbcseq$time1, pbcseq$time2, pbcseq$event)
py <- as.party(Rtree)
py$fitted[["(response)"]]<- survobj
plot(py)


