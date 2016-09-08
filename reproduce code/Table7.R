######################################################
## make sure to first load the necessary source code in the 
## "source code" file before running code below.
## If you are not sure which source code to load,
## load all of them just to be safe. 
######################################################

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

source("rpart_LTRC.R")
source("LTRC_ctree.R")
source("timevarying_type1.r")
source("timevarying_type3.r")
source("timevarying_type1_Ctnus.r")
source("timevarying_type3_Ctnus.r")

struct_recover <- function(n = 200, Dist = "Exponential", Jumptype = 1, Ctype = 1, Var = "Binary"){
	set.seed(1)

	Count_X1 = 0 
	Count_X2 = 0
	Count_S = 0
	
	Rount_X1 = 0 
	Rount_X2 = 0
	Rount_S = 0
	
	W=0
	while(W < 1000){
		if(Var == "Binary"){
			if(Jumptype == 1){
				Data1 <- try(Timevarying_gnrt(N = n, Distribution = Dist, censor.rate = Ctype), silent=TRUE)
				if(class(Data1) == "try-error"){
					next;
				}
			}else if(Jumptype == 3){
				Data1 <- try(Timevarying_gnrt3(N = n, Distribution = Dist, censor.rate = Ctype), silent=TRUE)
				if(class(Data1) == "try-error"){
					next;
				}
			}
		}else if(Var == "Continuous"){
			if(Jumptype == 1){
				Data1 <- try(Timevarying_gnrt_ctn(N = n, Distribution = Dist, censor.rate = Ctype), silent=TRUE)
				if(class(Data1) == "try-error"){
					next;
				}
			}else if(Jumptype == 3){
				Data1 <- try(Timevarying_gnrt3_ctn(N = n, Distribution = Dist, censor.rate = Ctype), silent=TRUE)
				if(class(Data1) == "try-error"){
					next;
				}
			}
		}else{stop("Wrong split variable type")}

		
		Tree.ctree <- LTRC_ctree(Surv(Start,Stop,Event)~X1+X2+X3+X4+X5, Data1)
		Tree.rpart <- rpart.LTRC(Surv(Start,Stop,Event)~X1+X2+X3+X4+X5, Data1)

		#check if LTRC ctree is correct 
		if( "X1" %in% getUsefulPredictors(Tree.ctree)){
			Count_X1 = Count_X1+1
		}
		if( "X2" %in% getUsefulPredictors(Tree.ctree)){
			Count_X2 = Count_X2+1
		}

		if( !(FALSE %in% (getPositions(Tree.ctree)== c("1","2","5"))) ){
			if( "X1" %in% getUsefulPredictors(Tree.ctree) && "X2" %in% getUsefulPredictors(Tree.ctree) ){
				if(length(Tree.ctree)==7 && width(Tree.ctree)==4 ){
					if(length(unique(getUsefulPredictors(Tree.ctree)))==2){
						Count_S = Count_S + 1;
					}
				}
			}
		}

		#check if LTRC rpart is correct 
		if( "X1" %in% unique(as.character(Tree.rpart$frame$var))){
			Rount_X1 = Rount_X1+1
		}
		if( "X2" %in% unique(as.character(Tree.rpart$frame$var))){
			Rount_X2 = Rount_X2+1
		}

		if(Tree.rpart$frame$var[1]=="X1"){
			if(Tree.rpart$frame$var[2]=="X2"){
				if(Tree.rpart$frame$var[5]=="X2"){
					if(length(unique(Tree.rpart$where))==4){
						Rount_S = Rount_S + 1;
					}
				}
			}
		}else if(Tree.rpart$frame$var[1]=="X2"){
			if(Tree.rpart$frame$var[2]=="X1"){
				if(Tree.rpart$frame$var[5]=="X1"){
					if(length(unique(Tree.rpart$where))==4){
						Rount_S = Rount_S + 1;
					}
				}
			}
		}
		W = W+1
	}##end of loop
	
	result <- list(Ctree = c(Count_X1, Count_X2, Count_S), Rtree = c(Rount_X1, Rount_X2, Rount_S))
	return(result)
}
###################=========================================================================
### n is the sample size, Dist specify the distribution
### Jumptype = 1 means dichotomous type I and Jumptype = 3 for type II
### Ctype = 0 means no censoring, Ctype = 1 means 20% censoring, and
###  Ctype = 0 means 50% censoring
###################=========================================================================
### Code below gives the results for sample size n=100 in type II case, i.e Table SM7
### Note that it gives the raw counts out of 1000 trials, divide by 10 
### to get the percentage rate report in paper
####################################################################################
Exp.0 <- struct_recover(n = 100, Dist = "Exponential", Jumptype = 3, Ctype = 0)
Exp.1 <- struct_recover(n = 100, Dist = "Exponential", Jumptype = 3, Ctype = 1)
Exp.2 <- struct_recover(n = 100, Dist = "Exponential", Jumptype = 3, Ctype = 2)

Web.0 <- struct_recover(n = 100, Dist = "Weibull", Jumptype = 3, Ctype = 0)
Web.1 <- struct_recover(n = 100, Dist = "Weibull", Jumptype = 3, Ctype = 1)
Web.2 <- struct_recover(n = 100, Dist = "Weibull", Jumptype = 3, Ctype = 2)

Gom.0 <- struct_recover(n = 100, Dist = "Gompertz", Jumptype = 3, Ctype = 0)
Gom.1 <- struct_recover(n = 100, Dist = "Gompertz", Jumptype = 3, Ctype = 1)
Gom.2 <- struct_recover(n = 100, Dist = "Gompertz", Jumptype = 3, Ctype = 2)

Exp.0
Exp.1
Exp.2
Web.0
Web.1
Web.2
Gom.0
Gom.1
Gom.2

