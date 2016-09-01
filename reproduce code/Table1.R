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
library(partykit)
library(rpart)
library(survival)
source("rpart_LTRC.R")
source("LTRC_ctree.R")
source("LTRC_structure_grt.r")

Struct_recover <- function(N = 200, distribution = "Exponential", Ctype = 1, Trunc = 2){
	set.seed(1)

	Count.ctree = 0 
	Count.rpart = 0
	
	for(W in 1:1000){
		DATA <- LTRC.generate(n = N, Dist = distribution, censor.type = Ctype, truncation = Trunc)

		Tree.ctree <- LTRC_ctree(Surv(Start,Obs,Event)~X1+X2+X3+X4+X5+X6, DATA)
		Tree.rpart <- rpart.LTRC(Surv(Start,Obs,Event)~X1+X2+X3+X4+X5+X6, DATA)

		#check if LTRC ctree is correct 
		if( !(FALSE %in% (getPositions(Tree.ctree)== c("1","2","5"))) ){
			if( sum(getUsefulPredictors(Tree.ctree) == c("X1","X2","X3"))==3 ){
				if(length(Tree.ctree)==7 && width(Tree.ctree)==4 ){
					Count.ctree = Count.ctree + 1;
				}
			}
		}

	      #check if LTRC rpart is correct 
		if(Tree.rpart$frame$var[1]=="X1"){
			if(Tree.rpart$frame$var[2]=="X2"){
				if(Tree.rpart$frame$var[5]=="X3"){
					if(length(unique(Tree.rpart$where))==4){
						Count.rpart = Count.rpart + 1;
					}
				}
			}
		}
	}##end of loop
	result <- c(Count.ctree, Count.rpart)
	return(result)
}

###################=========================================================================
########## Code below only gives the N=300 result,
########## Change N=100 and N=500 to get the complete result.
########## Note that it gives the row counts out of 1000 trials. Divide by 10 
########## to get the percentage as in the paper.
########## Ctype=1 means light censoring, while Ctype = 2 means heavy censoring
########## Trunc = 1 means the truncation time has U[0,1] distribution
########## Trunc = 3 means the truncation time has U[0,3] distribution

Exp.1.1 <- Struct_recover(N = 300, distribution = "Exponential", Ctype = 1, Trunc = 1)
Exp.2.1 <- Struct_recover(N = 300, distribution = "Exponential", Ctype = 2, Trunc = 1)
Exp.1.2 <- Struct_recover(N = 300, distribution = "Exponential", Ctype = 1, Trunc = 2)
Exp.2.2 <- Struct_recover(N = 300, distribution = "Exponential", Ctype = 2, Trunc = 2)
Exp.1.3 <- Struct_recover(N = 300, distribution = "Exponential", Ctype = 1, Trunc = 3)
Exp.2.3 <- Struct_recover(N = 300, distribution = "Exponential", Ctype = 2, Trunc = 3)

WebI.1.1 <- Struct_recover(N = 300, distribution = "Weibull-I", Ctype = 1, Trunc = 1)
WebI.2.1 <- Struct_recover(N = 300, distribution = "Weibull-I", Ctype = 2, Trunc = 1)
WebI.1.2 <- Struct_recover(N = 300, distribution = "Weibull-I", Ctype = 1, Trunc = 2)
WebI.2.2 <- Struct_recover(N = 300, distribution = "Weibull-I", Ctype = 2, Trunc = 2)
WebI.1.3 <- Struct_recover(N = 300, distribution = "Weibull-I", Ctype = 1, Trunc = 3)
WebI.2.3 <- Struct_recover(N = 300, distribution = "Weibull-I", Ctype = 2, Trunc = 3)


WebD.1.1 <- Struct_recover(N = 300, distribution = "Weibull-D", Ctype = 1, Trunc = 1)
WebD.2.1 <- Struct_recover(N = 300, distribution = "Weibull-D", Ctype = 2, Trunc = 1)
WebD.1.2 <- Struct_recover(N = 300, distribution = "Weibull-D", Ctype = 1, Trunc = 2)
WebD.2.2 <- Struct_recover(N = 300, distribution = "Weibull-D", Ctype = 2, Trunc = 2)
WebD.1.3 <- Struct_recover(N = 300, distribution = "Weibull-D", Ctype = 1, Trunc = 3)
WebD.2.3 <- Struct_recover(N = 300, distribution = "Weibull-D", Ctype = 2, Trunc = 3)

Log.1.1 <- Struct_recover(N = 300, distribution = "Lognormal", Ctype = 1, Trunc = 1)
Log.2.1 <- Struct_recover(N = 300, distribution = "Lognormal", Ctype = 2, Trunc = 1)
Log.1.2 <- Struct_recover(N = 300, distribution = "Lognormal", Ctype = 1, Trunc = 2)
Log.2.2 <- Struct_recover(N = 300, distribution = "Lognormal", Ctype = 2, Trunc = 2)
Log.1.3 <- Struct_recover(N = 300, distribution = "Lognormal", Ctype = 1, Trunc = 3)
Log.2.3 <- Struct_recover(N = 300, distribution = "Lognormal", Ctype = 2, Trunc = 3)

Bath.1.1 <- Struct_recover(N = 300, distribution = "Bathtub", Ctype = 1, Trunc = 1)
Bath.2.1 <- Struct_recover(N = 300, distribution = "Bathtub", Ctype = 2, Trunc = 1)
Bath.1.2 <- Struct_recover(N = 300, distribution = "Bathtub", Ctype = 1, Trunc = 2)
Bath.2.2 <- Struct_recover(N = 300, distribution = "Bathtub", Ctype = 2, Trunc = 2)
Bath.1.3 <- Struct_recover(N = 300, distribution = "Bathtub", Ctype = 1, Trunc = 3)
Bath.2.3 <- Struct_recover(N = 300, distribution = "Bathtub", Ctype = 2, Trunc = 3)
