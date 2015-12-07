# Least inverse squares M Estimator
mestimateMean<-function(MeansVector,ErrorsVector) {
	SquaredErrors<-ErrorsVector^2
	Numerator<-MeansVector/SquaredErrors
	InverseErrors<-sum(1/SquaredErrors)
	Mean<-sum(Numerator)/InverseErrors
	Error<-1/sqrt(InverseErrors)
	return(c(Mean,Error))
	}

# Dependency of normalOverlap()
findDensity<-function(x,Mean1,Mean2,StdDev1,StdDev2) {
    Density1<-dnorm(x,mean=Mean1,sd=StdDev1)
    Density2<-dnorm(x,mean=Mean2,sd=StdDev2)
    pmin(Density1,Density1)
	}

# Find the overlap of the mestimated distributions
normalOverlap<-function(MeansVector1,ErrorsVector1,MeansVector2,ErrorsVector2) {
	FirstDist<-mestimateMean(MeansVector1,ErrorsVector1)
	SecondDist<-mestimateMean(MeansVector2,ErrorsVector2)
	return(integrate(findDensity,-Inf,Inf,Mean1=FirstDist[1],Mean2=SecondDist[1],StdDev1=FirstDist[2],StdDev=FirstDist[3])
	}