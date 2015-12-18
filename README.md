# paleobiologyDatabase.R
R Functions for downloading, cleaning, culling, or analyzing fossil data from the Paleobiology Database. Developed and maintained by [Andrew Zaffos](https://macrostrat.org) as part of the [Paleobiology Database](https://paleobiodb.org) and [Macrostrat Database](https://macrostrat.org) tech development initiatives at the University of Wisconsin - Madison.

## Contents
+ [Creative Commons License](#creative-commons-license) # The creative commons license for modules in this repository.
+ [Version and Change Log](#version-and-change-log) # Information about changes to this repository.
+ [communityMatrix.R](#communitymatrixr) # Downloading and cleaning [Paleobiology Database](#https://paleobiodb.org) (PBDB) data, and making a community matrix.
+ [cullMatrix.R](#cullmatrixr) # Culling a communty matrix of depauperate samples and rare taxa.
+ [subsampleRichness.R](#subsamplerichnessr) # A set of subsampling functions for standardizing sampled taxonomic richness.
+ [partitionDiversity.R](#partitiondiversityr) # A set of functions for calculating alpha, beta, and gamma diversity of a community matrix.

## Creative Commons License
All code within the [paleobiologyDatabase.R](https://github.com/aazaff/paleobiologyDatabase.R) repository is covered under a Creative Commons License [(CC BY-NC 4.0)](http://creativecommons.org/licenses/by-nc/4.0/). This license requires attribution to [Andrew A. Zaffos](https://macrostrat.org/#contact) and [Steven M. Holland](https://strata.uga.edu), and does not allow for commercial use.

## Version and Change Log
This is v0.02 of the paleobiologyDatabase.R repository. The repository has four functional modules: [communityMatrix.R](#communitymatrixr), [cullMatrix.R](#cullmatrixr), [subsampleRichness.R](#subsamplerichnessr), and [partitionDiversity.R](#partitiondiversityr). Two incomplete modules are also currently uploaded: basicStatistics.R and gaussianOccupancy.R. These modules are still under development and their use is discouraged.

The next module will add support for several dual-concept diversity indices (e.g., True Shannon's Entropy) to some of these modules.

+ v.0.030 - Added [partitionDiversity.R](#partitiondiversityr) module.
+ v.0.025 - Upgraded [downloadPBDB( )](#downloadpbdb-) to use https instead of http.
+ v.0.024 - Removed support for [basicStatistics.R](#basicstatisticsr) module until additional functions come online.
+ v.0.023 - Removed  communityMatrix( ) and replaced it with the identical [presenceMatrix( )](#presencematrix-) function to make it more explicit that it is creating a presence-absence dataset. Added the [abundanceMatrix( )](#abundancematrix-) function, which makes a matrix with abundances.
+ v.0.022 - Added [basicStatistics.R](#basicstatisticsr) module. Currently only has one supported function, [mestimateMean( )](#mestimatemean-), which calculates the least inverse squares M-estimated mean and error. More functions coming soon.
+ v.0.021 - Added [resampleIndividuals( )](#resampleindividuals-) to [subsampleRichness.R](#subsamplerichnessr) module.
+ v.0.020 - Added [subsampleRichness.R](#subsamplerichnessr) module. Changed repository name from CleaningPBDB to paleobiologyDatabase.R. Added new function, [softCull( )](#softcull-), to cullMatrix.R module.
+ v.0.010 - Added [communityMatrix.R](#communitymatrixr) and [cullMatrix.R](#cullmatrixr) modules.

## communityMatrix.R
A set of functions for downloading data from the Paleobiology Database, and organizing it into a community matrix of samples by taxa, where "samples" are based on  one of the variables in the Paleobiology Database - e.g., Plate ID, Geologic Age.

Can be accessed directly in R using:

````
source("https://raw.githubusercontent.com/aazaff/paleobiologyDatabase.R/master/communityMatrix.R")
````

##### downloadPBDB( )
````
# Download data from PBDB by taxonomic group and geologic interval.

# Parameter Taxa must be a vector of one or more taxon names (as a character string), no default.
# Parameter StartInterval must be a single interval name accepted by the PBDB, default is "Pliocene"
# Parameter StopInterval must be a single interval name accepted by the PBDB, default is "Pleistocene" 

DataPBDB<-downloadPBDB(Taxa=c("Bivalvia","Gastropoda"),StartInterval="Cambrian",StopInterval="Pleistocene")
````

##### downloadTime( )
````
# Download Timescale definitions from Macrostrat.

# Parameter Timescale must be a timescale recognized by the macrostrat API, no default
# A list of Timescale defs can be seen here https://macrostrat.org/api/defs/timescales?all

Epochs<-downloadTime(Timescale="international epochs")
````

##### constrainAges( )
````
# Assign fossil occurrences to different ages, then remove occurrences that are not temporally 
# constrained to a single interval.

# Parameter DataPBDB is a dataset downloaded from the PBDB - i.e., using downloadPBDB( )
# Parameter Timescale is a dataset downloaded from Macrostrat - i.e., using downloadTime( )

ConstrainedAges<-constrainAges(DataPBDB=DataPBDB,Timescale=Epochs)
````

##### cleanGenus( )
````
# Cleans the genus field of the PBDB data by removing subgenera and NA's. This is an important step when
# working with genus level data in the PBDB, as the "genus" column often erroneously includes
# subgenus information.

# Parameter DataPBDB is a dataset downloaded from the PBDB - i.e., using downloadPBDB( )

CleanedPBDB<-cleanGenus(DataPBDB)
````

##### presenceMatrix( )
````
# Create a community matrix of samples v. species, using elements within one of the PBDB columns
# (e.g., geoplate, early_interval) as the definition of a sample. This creates a presence-absence
# matrix of 1's (presence) and 0's (absence).

# This is just a renamed version of the now deprecated function communityMatrix( ).

# Parameter DataPBBDB is a dataset downloaded from the PBDB - i.e., using downloadPBDB( )
# Parameter SampleDefinition is the column name defining samples

CommunityMatrix<-presenceMatrix(DataPBDB,SampleDefinition="geoplate")
````

##### abundanceMatrix( )
````
# Create a community matrix of samples v. species, using elements within one of the PBDB columns
# (e.g., geoplate, early_interval) as the definition of a sample. This creates an "abundance"
# matrix, which uses the number of occurrences a genus has within the "sample" as its abundance.
# Because the theoretical and operational meaning of occurrences in the Paleobiology Database is ill-defined
# I recommend using presenceMatrix( ) instead if possible.

# Parameter DataPBBDB is a dataset downloaded from the PBDB - i.e., using downloadPBDB( )
# Parameter SampleDefinition is the column name defining samples

CommunityMatrix<-abundanceMatrix(DataPBDB,SampleDefinition="geoplate")
````

## cullMatrix.R
A set of functions for removing depauperate and rare taxa from community matrices of samples by taxa.

Can be accessed directly in R using:

````
source("https://raw.githubusercontent.com/aazaff/paleobiologyDatabase.R/master/cullMatrix.R")
````
##### cullMatrix( )
````
# Cull a community matrix of depauperate samples and rare taxa. Written by S.M. Holland.

# Parameter x is a community matrix, no default.
# Parameter minOccurrences is the minimum number of occurrences for each taxon, default = 2
# Parameter minDiversity is the minimum number of taxa within each sample, default=2

CulledMatrix<-cullMatrix(CommunityMatrix,minOccurrences=5,minDiversity=5)
````

##### softCull( )
````
# A variant of cullMatrix( ) that returns NA when there are no samples or taxa left
# rather than throwing an error. Useful when culling multiple matrices in a loop, 
# and you would rather skip depauperate matrices than break the loop. 
# Not recommended otherwise.

# Parameter x is a community matrix, no default.
# Parameter minOccurrences is the minimum number of occurrences for each taxon, default = 2
# Parameter minDiversity is the minimum number of taxa within each sample, default=2

CulledMatrix<-softCull(CommunityMatrix,minOccurrences=5,minDiversity=5)
````

## subsampleRichness.R
Functions for standardizing taxonomic richness. The multicore versions use the [doMC](https://cran.r-project.org/web/packages/doMC/vignettes/gettingstartedMC.pdf) package and its dependencies.

Can be accessed directly in R using:

````
source("https://raw.githubusercontent.com/aazaff/paleobiologyDatabase.R/master/subsampleRichness.R")
````

##### subsampleEvenness( )
````
# A function that subsamples richness based on evenness. Often referred to as "Shareholder Quorum Subsampling".
# An optimized version of John Alroy's original function by Steven M. Holland.

# Parameter Abundance is a vector of abundances.
# Parameter Quota is a value between 0 and 1, the default is set to 0.9.
# Parameter Trials determines how many iterations of the bootstrap are performed, default = 100
# Parameter IgnoreSingletons determines whether or not to ignore singletons, default is FALSE.
# Parameter ExcludeDominant determines whether or not to ignore the most dominant taxon
# Excluding the abundant taxon is recommended by Alroy, but the default is set to FALSE.

SubsampledRichness<-subsampleEvenness(Abundance,Quota=0.5,Trials=100,IgnoreSingletons=FALSE,ExcludeDominant=FALSE)
````

##### multicoreEvenness( )
````
# A multicore version of subsampleEvenness( ). Be warned that multicoreEvenness( ) is not automatically
# faster than subsampleEvenness( ), particularly for a low number of trials. Its use is not recommended
# for small abundance vectors or a small numbers of trials. Requires the doMC package.

# Parameter Abundance is a vector of abundances.
# Parameter Quota is a value between 0 and 1, the default is set to 0.9.
# Parameter Trials determines how many iterations of the bootstrap are performed, default = 1000
# Parameter IgnoreSingletons determines whether or not to ignore singletons, default is FALSE.
# Parameter ExcludeDominant determines whether or not to ignore the most dominant taxon
# Excluding the abundant taxon is recommended by Alroy, but the default is set to FALSE.
# Parameter Cores sets the number of processor cores, default = 4.

SubsampledRichness<-multicoreEvenness(Abundance,Quota=0.5,Trials=100,IgnoreSingletons=FALSE,ExcludeDominant=FALSE,Cores=4)
````

##### subsampleIndividuals( )
````
# A function that subsamples richness based on a fixed number of individuals. Often referred to as "rarefaction".

# Parameter Abundance is a vector of abundances.
# Parameter Quota is the number of individuals to be subsampled. If the Quota is greater than
# the number of individuals, the function will print a warning and return the unstandardized richness.
# Parameter Trials determines how many iterations of the bootstrap are performed, default = 100

SubsampledRichness<-subsampleIndividuals(Abundance,Quota,Trials=100)
````

##### multicoreIndividuals( )
````
# A multicore version of subsampleIndividuals( ). Be warned that multicoreIndividuals( ) is not automatically
# faster than subsampleIndividuals( ), particularly for a low number of trials. Its use is not recommended
# for small abundance vectors or a small numbers of trials. Requires the doMC package.

# Parameter Abundance is a vector of abundances.
# Parameter Quota is the number of individuals to be subsampled. If the Quota is greater than
# the number of individuals, the function will print a warning and return the unstandardized richness.
# Parameter Trials determines how many iterations of the bootstrap are performed, default = 1000
# Parameter Cores sets the number of processor cores, default = 4.

SubsampledRichness<-multicoreIndividuals(Abundance,Quota,Trials=1000,Cores=4)
````

##### resampleIndividuals( )
````
# A specialized variant of subsampleIndividuals( ). If the quota is greater than the number of individuals
# it will switch to sampling with replacement. This allows for diversity in those samples to be lower
# than the quota. 

# Caution: This is a non-standard approach.

# Parameter Abundance is a vector of abundances.
# Parameter Quota is the number of individuals to be subsampled. 
# Parameter Trials determines how many iterations of the bootstrap are performed, default = 100

resampleIndividuals<-resampleIndividuals(Abundance,Quota,Trials=100)
````

## partitionDiversity.R
Functions for calculating alpha, beta, and gamma richness of a community matrix. See [communityMatrix.R](#communitymatrixr) for functions to make such a matrix with Paleobiology Database data and [cullMatrix.R](#cullmatrixr) for functions to cull and prepare
such a dataset. 

Some of these functions were presented in Holland, SM (2010) Additive diversity partitioning in palaeobiology: revisiting Sepkoski’s question. *Paleontology* 53:1237-1254. Namely, [taxonAlphaContributions( )](#taxonalphacontributions-), [taxonBetaContributions( )](#taxonbetacontributions-), and [sampleBetaContributions( )](#sampleBetaContributions-).

Other methods of beta calculation come from the equations presented in Tuomisto, H (2010) A diversity of beta diversities: straightening up a concept gone awry. Part 1. Defining beta diversity as a function of alpha and gamma diversity. *Ecography* 33:2-22.

Can be accessed directly in R using:

````
source("https://raw.githubusercontent.com/aazaff/paleobiologyDatabase.R/master/partitionDiversity.R")
````

##### taxonAlphaContributions( )
````
# Returns vector of each taxon’s contribution to alpha diversity. Written by S.M. Holland.

# Parameter x is a community matrix of presence-absence data.

TaxonAlpha<-taxonAlphaContributions(x=PresenceMatrix)
````

##### taxonBetaContributions( )
````
# Returns vector of each taxon’s contribution to beta diversity. 
# Be warned that if you are using a hierarchichal partitioning scheme # that this function *always* calculate between-sample beta,
# with sample being defined by your matrix. You must pre-aggregate samples # in the community matrix before you can calculate the beta # diversity of a higher level in the hierarchy. Written by S.M. Holland.

# Parameter x is a community matrix of presence-absence data.

TaxonBeta<-taxonBetaContributions(x=PresenceMatrix)
````

##### sampleBetaContributions( )
````
# Returns vector of each sample’s contribution to beta diversity. Written by S.M. Holland.

# Parameter x is a community matrix of presence-absence data.

TaxonBeta<-sampleBetaContributions(x=PresenceMatrix)
````

##### meanAlpha( )
````
# Returns mean alpha diversity (richness) of samples. Written by S.M. Holland.

# Parameter x is a community matrix of presence-absence data.

AlphaDiversity<-meanAlpha(x=PresenceMatrix)
````

##### beta( )
````
# Returns beta diversity (richness) of samples. Written by S.M. Holland.

# Parameter x is a community matrix of presence-absence data.

BetaDiversity<-beta(x=PresenceMatrix)
````

##### gamma( )
````
# Returns gamma (total) diversity (richness) of matrix. Written by S.M. Holland.

# Parameter x is a community matrix of presence-absence data.

GammaDiversity<-gamma(x=PresenceMatrix)
````

##### traditionalAlpha( )
````
# Calculate alpha diversity in the traditional manner, averaging sample richness.
# Should always be equal to meanAlpha( ) function.

# Parameter x is a community matrix of presence-absence data.

AlphaDiversity<-traditionalAlpha(x=PresenceMatrix)
````

##### traditionalBeta( )
````
# Calculate beta diversity in the traditional manner ADP manner Beta = Gamma - Alpha
# Should always be equal to beta( ) function.

# Parameter x is a community matrix of presence-absence data.

BetaDiversity<-traditionalBeta(x=PresenceMatrix)
````

##### multiplicativeBeta( )
````
# Calculate beta diversity in the traditional multiplicative manner. Beta = Gamma/Alpha

# Parameter x is a community matrix of presence-absence data.

BetaDiversity<-multiplicativeBeta(x=PresenceMatrix)
````

##### completeTurnovers( )
````
# Calculate Whittaker's effective species turnover, the number of complete effective species turnovers among samples in the dataset. 
# Beta = (Gamma-Alpha)/Alpha

# Parameter x is a community matrix of presence-absence data.

BetaDiversity<-completeTurnovers(x=PresenceMatrix)
````

##### proportionNonendemic( )
````
# Proportional effective species turnover, the proporition of species in the region not limited to a single sample - i.e., the 
# proportion of "non-endemic" taxa. Beta = (Gamma-Alpha)/Gamma

# Parameter x is a community matrix of presence-absence data.

BetaDiversity<-proportionNonendemic(x=PresenceMatrix)
````
