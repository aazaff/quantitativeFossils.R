# The Following are Steve Holland's Alpha, Beta, Gamma functions

# All functions assume that x is a presence-absence matrix, with
# columns corresponding to taxa and rows corresponding to sampling 
# units. The communityMatrix.R module on the paleobiologyDatabase.R
# repository provides functions for making such matrices.

# Data frame must be culled, such that all samples contain at 
# least one taxon, and all taxa occur in at least one sample.

# returns vector of each taxon’s contribution to alpha diversity
taxonAlphaContributions <- function(x) {
	nj <- apply(x, MARGIN=2, FUN=sum)
	N <- nrow(x)
	alphaj <- nj/N
	names(alphaj) <- colnames(x)
	return(alphaj)
	}

# returns vector of each taxon’s contribution to beta diversity
taxonBetaContributions <- function(x) {
	nj <- apply(x, MARGIN=2, FUN=sum)
	N <- nrow(x)
	nminusj <- N - nj
	betaj <- nminusj / N
	names(betaj) <- colnames(x)
	return(betaj)
	}

# returns vector of each sample’s contribution to beta diversity
sampleBetaContributions <- function(x) {
	betaj <- taxonBetaContributions(x)
	nj <- apply(x, MARGIN=2, FUN=sum)
	numSamplingUnits <- nrow(x)
	betai <- vector(mode="numeric", length=numSamplingUnits)
	for (i in 1:numSamplingUnits) {
		betai[i] <- sum(betaj / nj * x[i, ])
		}
	names(betai) <- rownames(x)
	return(betai)
	}

# returns mean alpha diversity of samples
meanAlpha <- function (x) {
	return(sum(taxonAlphaContributions(x)))
	}

# returns beta diversity among samples
beta <- function (x) {
	return(sum(taxonBetaContributions(x)))
	}

# returns gamma diversity of all samples combined
gamma <- function(x) {
	return(ncol(x))
	}

# Calculate alpha diversity in the traditional manner, average sample richness
traditionalAlpha<-function(x) {
	return(mean(apply(x,1,sum)))
	}

# Calculate beta diversity in the traditional manner, gamma - alpha
traditionalBeta<-function(x) {
	Alpha<-traditionalAlpha(x)
	Gamma<-gamma(x)
	Beta<-Gamma-Alpha
	return(Beta)
	}

# "True local diversity ratio" of Tuomisto 2010
# This quantifies how many times as rich in effective species gamma is than alpha
multiplicativeBeta<-function(x) {
	Alpha<-traditionalAlpha(x)
	Gamma<-gamma(x)
	Beta<-Gamma/Alpha
	return(Beta)
	}

# Whittaker's effective species turnover, the number of complete effective species
# turnovers among compositional units in the dataset
completeTurnovers<-function(x) {
	Alpha<-traditionalAlpha(x)
	Gamma<-gamma(x)
	Beta<-(Gamma-Alpha)/Alpha
	return(Beta)
	}

# Proportional effective species turnover, the proporition of species in the region
# not limited to a single sample - i.e., the proportion of non-endemic taxa.
proportionNonendemic<-function(x) {
	Alpha<-traditionalAlpha(x)
	Gamma<-gamma(x)
	Beta<-(Gamma-Alpha)/Gamma
	return(Beta)
	}