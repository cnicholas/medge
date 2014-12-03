#' Decompose the difference between two weighted averages.  This method will calculate the weighted averages 
#' for two time periods and execute the Marshall Edgeworth method. 
#' 
#' @param groups A vector of labels for the groups being compared.
#' @param numerator.t1 A vectors of values for the numerator of the calculation for t1.
#' @param numerator.t2 A vectors of values for the numerator of the calculation for t2.
#' @param denominator.t1 A vectors of values for the denominator of the calculation for t1.
#' @param denominator.t2 A vectors of values for the denominator of the calculation for t2.
#' @param scale.rates.by value to scale rates by (default=1 though 100 improves legibility of effects) 
#' @return object of class medge
#' @export
#' @examples
#' medge(
#' 			groups=c("A", "B", "C", "D", "E")
#' 			,numerator.t1=c(10,20,30,40,50)
#' 			,numerator.t2=c(11,19,35,35,52)
#' 			,denominator.t1=c(20,40,60,80,100)
#' 			,denominator.t2=c(30,50,50,90,90)
#' 			)
medge <- function(groups,numerator.t1,numerator.t2,denominator.t1,denominator.t2, scale.rates.by=1){
	
	#Validations
	#Length of all args must be equal
	equal.arg.lengths = all(c(length(numerator.t1), length(numerator.t2), length(denominator.t1), length(denominator.t2))==length(groups))
	if (!equal.arg.lengths) stop("Argument lengths must be equal")
	
	#Validate Divide by Zero non-zero numerator / zero denonminator
	rates.t1 <- (numerator.t1 / denominator.t1)
	rates.t2 <- (numerator.t2 / denominator.t2)
	
	if(any(is.infinite(rates.t1))) stop("Divide by zero ocurred during t1 rate calculation.",call.=FALSE, domain=NA)
	if(any(is.infinite(rates.t2))) stop("Divide by zero ocurred during t2 rate calculation.",call.=FALSE, domain=NA)
	
	#Update 0/0 situations to result in rate of 0 (Occurs when a new group enters or exits the population)
	rates.t1 <- ifelse(is.nan(rates.t1), 0, rates.t1)
	rates.t2 <- ifelse(is.nan(rates.t2), 0, rates.t2)
	
	#Create medge object
	object <- list(
			  groups=groups
			, numerators=list(t1=numerator.t1, t2=numerator.t2)
			, denominators=list(t1=denominator.t1,t2=denominator.t2)
	 		, rates=list(
					t1 = rates.t1 * scale.rates.by
					, t2 = rates.t2 * scale.rates.by
			)
			, overall.rates = list(
					t1 = sum(numerator.t1)/sum(denominator.t1) * scale.rates.by
					, t2 = sum(numerator.t2)/sum(denominator.t2) * scale.rates.by
			)
			, overall.weights=list(
					t1 = sum(denominator.t1)
					,t2 = sum(denominator.t2)
			)
	)
	
	object$weights=list(
			t1 = (denominator.t1/object$overall.weights$t1)
			, t2 = (denominator.t2/object$overall.weight$t2)
	)
	object$rate.effects = (object$rates$t2-object$rates$t1)*((object$weights$t2+object$weights$t1)/2)
	object$weight.effects = (((object$rates$t2-object$overall.rates$t2)+(object$rates$t1-object$overall.rates$t1))/2) * (object$weights$t2-object$weights$t1) 
	object$rankings=rank(abs(object$rate.effects+object$weight.effects))
	
	object$overall.rate.effect = sum(object$rate.effects)
	object$overall.weight.effect = sum(object$weight.effects)
	object$global.effect = object$overall.rate.effect + object$overall.weight.effect
	
	class(object) <- "medge"
	return(object)
}

#' Decompose the difference between two weighted averages.  This method will calculate the weighted averages for two time periods and execute the Marshall Edgeworth method. 
#' 
#' @param x medge
#' @param ... additional parameters
#' @return str(x)
#' @export
print.medge <- function(x, ...) str(x)

#' Provides a summary of the differences between the weighted averages 
#' 
#' @param object medge
#' @param digits digits
#' @param ... additional parameters
#' @export
summary.medge <- function(object, digits=getOption["digits"], ...){
	cat("The difference in rates between t1 (", object$overall.rates$t1, ") and t2 (", object$overall.rates$t2, ") is: ", object$overall.rates$t2 - object$overall.rates$t1,"\r\n",sep="")
	cat("The rate effect accounts for ", object$overall.rate.effect, "percentage points of the ",object$overall.rates$t2 - object$overall.rates$t1," difference.\r\n")
	cat("The weight effect accounts for ", object$overall.weight.effect, "percentage points of the change",object$overall.rates$t2 - object$overall.rates$t1," difference.\r\n")
}

#' Convert the medgworth object into a data.frame 
#' 
#' @param object medge
#' @param sort.by.ranking if true sort by the ranking of the global effect
#' @return data.frame of the marshall edgeworth decomposition 
#' @export
as.data.frame.medge <- function(object, sort.by.ranking=FALSE,...){
	
	opt_scipen <- options("scipen")
	options(scipen=999)
	
	retval=data.frame(
			groups=object$groups
			,numerator.t1=object$numerators$t1
			,numerator.t2=object$numerators$t2
			,denominator.t1=object$denominators$t1
			,denominator.t2=object$denominators$t2
			,rates.t1 = object$rates$t1
			,rates.t2 = object$rates$t2
			,weights.t1 = object$weights$t1
			,weights.t2 = object$weights$t2
			,ranking=object$rankings
			,weight.effects=object$weight.effects  
			,rate.effects=object$rate.effects
	)
	if(sort.by.ranking){
		retval=retval[order(-retval$ranking),]
	}
	print(head(retval))
	print(".....")
	options(scipen=opt_scipen[[1]])
	retval
}
#' Plot a bar graph of the overall rate and weight effects 
#' 
#' @param object
#' @param ... additional parameters
#' @export
plot.medge<-function(object,...){
	barplot(c(object$overall.rate.effect,object$overall.weight.effect), main=paste("Marshall Edgeworth Decomposition\r\n","Explanation for the: ",round(object$global.effect,3), " change.",sep=""), 
			xlab="Contributions",names.arg=c("Rate Effect","Weight Effect"),ylab="Percentage Points")
}

