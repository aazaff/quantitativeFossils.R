# Least inverse squares M Estimator
mestimateMean<-function(MeansVector,ErrorsVector) {
	SquaredErrors<-ErrorsVector^2
	Numerator<-MeansVector/SquaredErrors
	InverseErrors<-sum(1/SquaredErrors)
	Mean<-sum(Numerator)/InverseErrors
	Error<-1/sqrt(InverseErrors)
	return(c(Mean,Error,length(MeansVector)))
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
	return(integrate(findDensity,-Inf,Inf,Mean1=FirstDist[1],Mean2=SecondDist[1],StdDev1=FirstDist[2],StdDev=SecondDist[2]))
	}

# Conduct a t test using the m-estimated means and errors
mestimateTtest<-function(MeansVector1,ErrorsVector1,MeansVector2,ErrorsVector2) {
	FirstDist<-mestimateMean(MeansVector1,ErrorsVector1)
	SecondDist<-mestimateMean(MeansVector2,ErrorsVector2)
	Numerator<-FirstDist[1]-SecondDist[2]
	DenomRight<-(1/FirstDist[3]+1/SecondDist[3])
	DenomTop1<-(FirstDist[3]-1)*FirstDist[2]^2
	DenomTop2<-(SecondDist[3]-1)*SecondDist[2]^2
	DenomTop<-DenomTop1+DenomTop2
	DenomBottom<-FirstDist[3]+SecondDist[3]-2
	DenomLeft<-DenomTop/DenomBottom
	Denominator<-sqrt(DenomLeft*DenomRight)
	t<-Numerator/Denominator
	}