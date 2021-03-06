library(rpart)
library(survival)

rpart.LTRC <- function(formula, data, weights=NULL, subset = NULL, no.SE = 0, control = rpart.control(cp = 0.001), cost){

	if(missing(data)){
		y <- eval(formula[[2]])
		predictors <- eval(formula[[3]])

		if (!inherits(y, "Surv") | length(as.list(formula[[2]]))!=4){
			stop(" Response must be a 'survival' object with 'Surv(time1, time2, event)' ")
		}

		Status <- y[,3L]
		Times <- y[,2L]

		##unique death times
		unique.times <- sort(unique(Times[Status == 1])) 
		
		temp <- coxph(y ~ 1)
		cumhaz.table <- basehaz(temp)
		
		#Check if Inf hazard exists
		if(sum(is.infinite(cumhaz.table$hazard))!=0){
			cumhaz.table <- subset(cumhaz.table, hazard != Inf)
		}

		cumhaz.table2 <- subset(cumhaz.table, time %in% unique.times) 
		
		cumhaz.times <- c(0, cumhaz.table2$time[-length(cumhaz.table2$time)], max(Times))
		cumhaz <- c(0, cumhaz.table2$hazard)
	
		Start.cumhaz <- approx(cumhaz.times, cumhaz, y[,1L])$y
		End.cumhaz <- approx(cumhaz.times, cumhaz, y[,2L])$y
	
		Newtime <- End.cumhaz - Start.cumhaz
		Formula = formula(paste(c( paste("cbind(Newtime,","Status)",sep = ""), "predictors"), collapse = "~"))
		result <- rpart(formula = Formula, method = "poisson", weights = weights, subset = subset, control = control, cost = cost)
		return(result)
	}else{
		Data <- data
		###if in the form of Surv(time1,time2,event)~predictors with data following
   		Response <- formula[[2]]
		predictors <- formula[[3]]
	
		if(length(as.list(Response))!=4){
			stop(" Response must be a 'survival' object with 'Surv(time1, time2, event)' ")
		}
		y.names <- c(as.character(as.list(Response)[[2]]),as.character(as.list(Response)[[3]]),as.character(as.list(Response)[[4]])) 
     	y.IDs <- match(y.names, names(Data))
		y <- Surv(Data[,y.IDs[1]],Data[,y.IDs[2]],Data[,y.IDs[3]])

		Status <- y[,3L]
		Times <- y[,2L]

		##unique death times
		unique.times <- sort(unique(Times[Status == 1])) 
		temp <- coxph(y ~ 1)
		cumhaz.table <- basehaz(temp)
		
		#Check if Inf hazard exists
		if(sum(is.infinite(cumhaz.table$hazard))!=0){
			cumhaz.table <- subset(cumhaz.table, hazard != Inf)
		}
		
		cumhaz.table2 <- subset(cumhaz.table, time %in% unique.times) 

		cumhaz.times <- c(0, cumhaz.table2$time[-length(cumhaz.table2$time)], max(Times))
		cumhaz <- c(0, cumhaz.table2$hazard)
	
		Start.cumhaz <- approx(cumhaz.times, cumhaz, y[,1L])$y
		End.cumhaz <- approx(cumhaz.times, cumhaz, y[,2L])$y

		Data$Newtime <- End.cumhaz - Start.cumhaz
      	Formula = formula(paste(c( paste("cbind(Newtime,",y.names[3],")",sep = ""), predictors), collapse = "~"))
  		DATA = Data[,-y.IDs[1:2]]
		Fit.tree <- rpart(formula = Formula, data = DATA, method = "poisson", weights = weights, subset = subset, control = control, cost = cost)
		
		##pruning the tree
		if( length(unique(Fit.tree$where))!=1 ){
			cventry <- which.min(Fit.tree$cptable[, "xerror"])

			if(no.SE == 0){
				cpcv <- Fit.tree$cptable[cventry, "CP"]
				Fit.tree <- prune(Fit.tree,cp=cpcv)
			}else{
				xerrorcv <- Fit.tree$cptable[cventry, "xerror"]
                 		sexerrorcv <- xerrorcv + Fit.tree$cptable[cventry,"xstd"] * no.SE
                  	cpcvse <- Fit.tree$cptable[which.max(Fit.tree$cptable[, "xerror"] <= sexerrorcv), "CP"]
                  	Fit.tree <- prune(Fit.tree, cp = cpcvse)
			}
		}
		return(Fit.tree)
	}
}

 
