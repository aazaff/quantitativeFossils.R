# From Marshall's (1990) adaptation of Strauss and Sadler.
# Brackets the extinction date
extinctionConfidenceSS <- function(OccurrenceAges, ConfidenceLevel=0.95)  {
    # Find the number of unique "Horizons"
    NumOccurrences<-length(unique(OccurrenceAges))-1
    Alpha<-((1-ConfidenceLevel)^(-1/NumOccurrences))-1
    Lower<-min(OccurrenceAges)
    Upper<-min(OccurrenceAges)-(Alpha*10)
    return(setNames(c(Lower,Upper),c("Earliest","Latest")))
    }

# A dependency of iterateParametersABM
integrand.neglambdas<-function(Lambda,Theta,OccurrenceAges,PRMean=0,PRStdDev=2)  {
    FinalVector<-vector("numeric",length(Lambda))
    for(i in 1:length(Lambda)) { 
        FinalVector[i]<-(sum(log((1-Lambda[i])/Theta * 1/(1-OccurrenceAges/Theta)^Lambda[i])))
        }
    FinalVector<-FinalVector+FinalVector(Lambda,PRMean,PRStdDev,log=TRUE) + log(1/Theta)
    FinalVector<-exp(FinalVector)
    return(output)
    }

  integrand.poslambdas <- function(L,th,x)  {
  # L = vector of lambdas,  th = theta value (scalar),  x = vector of strat. positions
    k <- length(L)    
    output <- rep(NA, k)
    for(i in 1:k) 
      output[i] <- ( sum( log( (1+L[i])/th * (x/th)^L[i] ) ) )
    output <- output + dnorm(L, prmean,prSD, log=T) + log(1/th)
    output <- exp(output)
    return(output)
  }

  integrand.thetasnegL <- function(th,L,x)  {
  # th = vector of thetas,  L = lambda value (scalar), x = vector of strat. positions
    k <- length(th) 
    output <- rep(NA, k)
    for(i in 1:k) 
      output[i] <- ( sum( log( (1-L)/th[i] * 1/(1-x/th[i])^L ) ) )
    output <- output + dnorm(L, prmean,prSD, log=T) + log(1/th)
    output <- exp(output)
    return(output)  
  }

  integrand.thetasposL <- function(th,L,x)  {
  # th = vector of thetas,  L = lambda value (scalar), x = vector of strat. positions
    k <- length(th) 
    output <- rep(NA, k)
    for(i in 1:k) 
      output[i] <- ( sum( log( (1+L)/th[i] * (x/th[i])^L) ) )
    output <- output + dnorm(L, prmean,prSD, log=T) + log(1/th)
    output <- exp(output)
    return(output)  
  }

  
  # the function below is used in likelihood calculations
  drefbeta <- function(x,L)  {
    if(L<=0)  { 
      return(dbeta(x, 1, 1-L))
    }  else  return(dbeta(x, 1+L, 1))
  }

# A dependency of extinctionConfidenceABM
iterateParametersABM<-function(OccurrenceAges,UpperTheta=500,LengthTheta=1000,LowerLambda=-10,UpperLamda=10,LengthLambda=40) {
    LambdaValues<-seq(LowerLamda,UpperLamda,length.out=LengthLamda)
    LambdaDensity<-vector("numeric",length(Lambda))
    ThetaValues<-seq(max(OccurrenceAges),UpperTheta,length.out=LengthTheta)
    ThetaDensity<-vector("numeric",length(LengthTheta))
    
    # increment lambda values, integrating over theta values for each
    for(i in 1:LengthLambda)  {
        LambdaDensity[i]<-ifelse(LambdaValues[i]<=0, 
            integrate(integrand.thetasnegL, max(OccurrenceAges),Theta=UpperTheta, Lambda=LambdaValues[i], OccurrenceAges=OccurrenceAges)$value,
            integrate(integrand.thetasposL, max(OccurrenceAges),Theta=UpperTheta, Lambda=LambdaValues[i], OccurrenceAges=OccurrenceAges)$value
            )
        }
        
        
        
        
        
        
    LambdaDensity<-LambdaDensity/sum(LambdaDensity) # normalize lambda pdf to unit area  
    # calculate posterior quantities
    # cutoff <- which.max( cumsum(Ldens) >= .5 )       # posterior median
    # Lmed <- Lvals[cutoff]
    Lmean <- sum(Lvals*Ldens)
    Lhat <- Lmean


# A simplified version of Steve Wang's adaptive beta method of extinction estimation
# See his supplementary files on DataDryad for a more comprehensive version with more options
extinctionConfidenceABM <- function(OccurrenceAges, ConfidenceLevel=0.95)    {
    # Find the starting age
    Base<-min(OccurenceAges)
    OccurenceAges<-OccurrenceAges-Base
    OccurrenceAges<-OccurrenceAges[-which.min(OccurrenceAges)]  # Remove the Base age from the pool of ages
    
    # Scale data so theta is approx. 100 (solely for numerical stability)
    MaxAge<-max(OccurrenceAges)
    NumOccurrences<-length(OccurrenceAges)
    Theta<-(NumOccurrences+1)/NumOccurrences*MaxAge
    ScaleFactor<-100/Theta
    OccurrenceAges<-OccurrenceAges*ScaleFactor
    MaxAge<-max(OccurrenceAges)
    

    
    
