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
	alphaj
	}

# returns vector of each taxon’s contribution to beta diversity
taxonBetaContributions <- function(x) {
	nj <- apply(x, MARGIN=2, FUN=sum)
	N <- nrow(x)
	nminusj <- N - nj
	betaj <- nminusj / N
	names(betaj) <- colnames(x)
	betaj
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
	betai
	}

# returns mean alpha diversity of samples
meanAlpha <- function (x) {
	sum(taxonAlphaContributions(x))
	}

# returns beta diversity among samples
beta <- function (x) {
	sum(taxonBetaContributions(x))
	}

# returns gamma diversity of all samples combined
gamma <- function(x) {
	ncol(x)
	}