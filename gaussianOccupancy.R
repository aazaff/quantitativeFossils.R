# Steve Holland Function for caluclating GLR
tBLparamsForSpeciesMatrix <- function(sampleValues, abundanceMatrix) {
	numTaxa <- ncol(abundanceMatrix)
	opt <- vector(mode="numeric", length=numTaxa)
	tol <- vector(mode="numeric", length=numTaxa)
	pmax <- vector(mode="numeric", length=numTaxa)
	for (t in 1:numTaxa) {
		x <- sampleValues[,1]
		y <- abundanceMatrix[,t]
		y[y>0] <- 1
		logitReg <- glm(y ~ x + I(x^2), family=binomial)
		b0 <- logitReg$coefficients[1]
		b1 <- logitReg$coefficients[2]
		b2 <- logitReg$coefficients[3]
		opt[t] <- (-b1)/(2*b2)
		tol[t] <- 1 / sqrt(-2*b2)
		pmax[t] <- 1 / (1 + exp(b1^2 / (4 * b2) - b0)) 
		}
	tBLParams <- data.frame(opt, tol, pmax)
	tBLParams$opt[is.na(tBLParams$tol)] <- NaN
	tBLParams$pmax[is.na(tBLParams$tol)] <- NaN
	return(tBLParams)
	}

# Calculate GLR Parameters
paramsGLR<-function (Abundances) {
	Latitudes<-as.numeric(rownames(Abundances))
	Latitudes<-cbind(Latitudes,c(1:length(Latitudes)))
	ParamsGLR<-tBLparamsForSpeciesMatrix(Latitudes,Abundances)
	rownames(ParamsGLR)<-colnames(Abundances)
	ParamsGLR<-na.omit(ParamsGLR)
	return(ParamsGLR)
	}

cullParams<-function (ParamsGLR,MinLat=-90,MaxLat=90,MaxRange=181) {
	for (i in 1:nrow(ParamsGLR)) {
		if (ParamsGLR[i,"pmax"]<0.001) {
			ParamsGLR[i,]<-NA
			}
		else if (ParamsGLR[i,"pmax"]>0.999) {
			ParamsGLR[i,]<-NA
			}
		else if (ParamsGLR[i,"opt"]<=MinLat) {
			ParamsGLR[i,]<-NA
			}
		else if (ParamsGLR[i,"opt"]>=MaxLat) {
			ParamsGLR[i,]<-NA
			}
		else if (ParamsGLR[i,"tol"]>=MaxRange) {
			ParamsGLR[i,]<-NA
			}
		}
	return(ParamsGLR)
	}
	
# Convert preferences into a quadratic formula
createQuadratic<-function(CreatedPreferences) {
	ParamMatrix<-matrix(NA,nrow=nrow(CreatedPreferences),ncol=3)
	colnames(ParamMatrix)<-c("b2","b1","b0")
	rownames(ParamMatrix)<-rownames(CreatedPreferences)
	ParamMatrix[,"b2"]<-(1/(CreatedPreferences[,2]^2))/-2
	ParamMatrix[,"b1"]<-(-(2*CreatedPreferences[,1]*ParamMatrix[,1]))
	LeftTerm<-log(CreatedPreferences[,3])
	MidTerm<-ParamMatrix[,2]*CreatedPreferences[,1]
	ThirdTerm<-ParamMatrix[,1]*(CreatedPreferences[,1]^2)
	ParamMatrix[,"b0"]<-LeftTerm-MidTerm-ThirdTerm
	return(ParamMatrix)
	}

# Find the abundance
abundanceLatitude<-function(QuadParams,Lower=-90,Upper=90,Breaks=1) {
	Lats<-seq(Lower,Upper,Breaks)
	FinalMatrix<-matrix(NA,nrow=nrow(QuadParams),ncol=length(Lats))
	colnames(FinalMatrix)<-Lats
	rownames(FinalMatrix)<-rownames(QuadParams)
	for (i in 1:nrow(FinalMatrix)) {
		Abundances<-(QuadParams[i,"b0"])+(QuadParams[i,"b1"]*Lats)+(QuadParams[i,"b2"]*Lats^2)
		FinalMatrix[i,]<-exp(Abundances)
		}	
	return(FinalMatrix)
	}
