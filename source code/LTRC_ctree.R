###require packages: partykit, survival###
library(partykit)
library(survival)

LTRC_ctree <- function(Formula, Data, Control = ctree_control()){
	Response <- as.character(Formula)[[2]]

	## logrank transformation for left-truncated and right-censored data
	## x2 is Surv(Left,right,event) object
	.logrank_trafo <- function(x2){

		unique.times <- unique(x2[,2][which(x2[,3]==1)])

		##event number and risk set at each unique times
		D <- rep(NA,length(unique.times))
		R <- rep(NA,length(unique.times))

		for(j in 1:length(unique.times)){
			D[j] = sum(unique.times[j]==x2[,2])
		}

		for(k in 1:length(unique.times) ){
 			value <- unique.times[k]
			R[k]<- sum(apply(x2[,1:2], 1, function(interval){interval[1] < value & value <= interval[2]}))
		}

		Ratio <- D/R

		Ratio <- Ratio[order(unique.times)]
		Nelson.Aalen <- cumsum(Ratio)
		Event.time <- unique.times[order(unique.times)]
	
		##cumulative hazard at left truncation time
		Left <- sapply(x2[,1],function(t){if(t< min(Event.time)) return(0) else return(Nelson.Aalen[max(which(Event.time <= t))])})
		Right <- sapply(x2[,2],function(t){if(t< min(Event.time)) return(0) else return(Nelson.Aalen[max(which(Event.time <= t))])})

		result<- x2[,3]-(Right-Left)

		return(matrix(as.double(result), ncol = 1))
	}
	##--------------------------------------------------------------
	h2 <- function(data, weights) {
		s <- data[, Response]
		s <- .logrank_trafo(s[weights > 0,])
		r <- rep(0, nrow(data))
		r[weights > 0] <- s
		matrix(as.double(r), ncol = 1)
	}
	#--------------------------------------------------------------------
	partykit::ctree(formula = Formula, data=Data, ytrafo=h2, control = Control)
}

########### An example, where the data name is DATA, Start is left truncation time, 
########### Obs is right censoring/observed time, Event is indicator with 1 or 0
###########---------- Tree <- LTRC_ctree(Surv(Start, Obs, Event)~X1+X2+X3, Data = DATA) -------------

