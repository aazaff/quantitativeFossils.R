# Steve Holland's Culling Functions
cullMatrix <- function(x, minOccurrences=2, minDiversity=2) {
	y <- x
	while (diversityFlag(y, minDiversity) | occurrencesFlag(y, minOccurrences)) {
		y <- cullTaxa(y, minOccurrences)
		y <- cullSamples(y, minDiversity)
	}
	y
}

# Dependency of cullMatrix()
cullTaxa <- function(x, minOccurrences) {
	PA <- x
 	PA[PA>0] <- 1
 	occurrences <- apply(PA, MARGIN=2, FUN=sum)
 	aboveMinimumOccurrences <- occurrences >= minOccurrences
 	y <- x[ ,aboveMinimumOccurrences]
 	if (length(y)==0) {print("no taxa left!")}
 	y
}

# Dependency of cullMatrix()
cullSamples <- function(x, minDiversity) {
	PA <- x
 	PA[PA>0] <- 1
 	diversity <- apply(PA, MARGIN=1, FUN=sum)
 	aboveMinimumDiversity <- diversity >= minDiversity
 	y <- x[aboveMinimumDiversity, ]
 	if (length(y[,1])==0) {print("no samples left!")}
 	y
}

# Dependency of cullMatrix()
occurrencesFlag <- function(x, minOccurrences) {
  	PA <- x
 	PA[PA>0] <- 1
 	occurrences <- apply(PA, MARGIN=2, FUN=sum)
 	if (min(occurrences) < minOccurrences) {
 		flag <- 1
 	}
 	else {
	 	flag <- 0
 	}
 	flag
}

# Dependency of cullMatrix()
diversityFlag <- function(x, minDiversity) {
 	PA <- x
 	PA[PA>0] <- 1
 	diversity <- apply(PA, MARGIN=1, FUN=sum)
 	if (min(diversity) < minDiversity) {
 		flag <- 1
 	}
 	else {
 		flag <- 0
 	}
 	flag
}

# Alternative to cullMatrix( ) that does not thrown an error, but returns a single NA
softCull <- function(x, minOccurrences=2, minDiversity=2) {
	y <- x
	while (diversityFlag(y, minDiversity) | occurrencesFlag(y, minOccurrences)) {
		y <- softTaxa(y, minOccurrences)
		if (is.na(y)==TRUE) {
			return(NA)
			}
		y <- softSamples(y, minDiversity)
		if (is.na(y)==TRUE) {
			return(NA)
			}
		}
	return(y)
	}
	
# Dependency of softCull()
softTaxa <- function(x, minOccurrences) {
	PA <- x
 	PA[PA>0] <- 1
 	occurrences <- apply(PA, MARGIN=2, FUN=sum)
 	aboveMinimumOccurrences <- occurrences >= minOccurrences
 	y <- x[ ,aboveMinimumOccurrences]
 	if (length(y)==0) {
 		y<-NA;
 		return(NA)
 		}
 	return(y)
	}

# Dependency of softCull()
softSamples <- function(x, minDiversity) {
	PA <- x
 	PA[PA>0] <- 1
 	diversity <- apply(PA, MARGIN=1, FUN=sum)
 	aboveMinimumDiversity <- diversity >= minDiversity
 	y <- x[aboveMinimumDiversity, ]
 	if (length(y[,1])==0) {
 		y<-NA;
 		return(NA)
 		}
 	return(y)
	}
