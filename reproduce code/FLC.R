## make sure to set the directory correctly 
setwd("C:/Users/wfu/Desktop/biostatistics/code")

library(survival)
library(rpart.plot)
source("rpart_LTRC.R")
source("LTRC_ctree.R")


#Adjust data & clean data
Data <- flchain
Data <- Data[!is.na(Data$creatinine),]
Data$End <- Data$age+Data$futime/365
DATA <- Data[Data$End > Data$age,]
names(DATA)[6] <- "FLC"

set.seed(1)
## Fit LTRCIT and regular survival trees
LTRC_Tree <-  LTRC_ctree(Surv(age, End, death)~sex+FLC+creatinine,DATA)
Ctree <- partykit::ctree(Surv(End, death)~sex+FLC+creatinine,DATA)


##Fit LTRCART and regular rpart trees
Tree.rpart <- rpart.LTRC(Surv(age, End, death)~sex+FLC+creatinine,DATA)
rpart <- rpart(Surv(End, death)~sex+FLC+creatinine,DATA)
## prune the rpart tree 
cventry <- which.min(rpart$cptable[, "xerror"])
cpcv <- rpart$cptable[cventry, "CP"]
rpart <- prune(rpart,cp=cpcv)

##############======Old Plots==========#############################
##### Old plots. Fig.9, Fig.10 in the first manuscript---############
plot(LTRC_Tree)
plot(Ctree)
rpart.plot.version1(Tree.rpart)
rpart.plot.version1(rpart)

##### Old Cox model results. Table 2, Table 3 in the first manuscript---############
##### Fit Cox model and obtain its results
LTRC_Cox<- coxph(Surv(age, End, death)~sex+FLC+creatinine,DATA)
Strat_Cox <- coxph(Surv(age, End, death)~sex:FLC+creatinine,DATA)

summary(LTRC_Cox)
summary(Strat_Cox)
#####################################################################



#####======New plots in the revised papaer===========#########

pdf("LTRC_FLC.pdf",width=10)
plot(LTRC_Tree)
dev.off()

pdf("Ctree_FLC.pdf",width=10)
plot(Ctree)
dev.off()

#####Convert the LTRCART object to compatible with ctree plot###############
survobj <- Surv(DATA$age, DATA$End, DATA$death)
py <- as.party(Tree.rpart)
py$fitted[["(response)"]]<- survobj

pdf("FLC_LTRC_rpart.pdf",width=10)
plot(py)
dev.off()

pdf("FLC_rpart.pdf",width=10)
plot(as.party(rpart))
dev.off()
