# paleobiologyDatabase.R
R Functions for downloading, cleaning, culling, or analyzing fossil data from the Paleobiology Database. Developed and maintained by [Andrew Zaffos](macrostrat.org) as part of the [Paleobiology Database](paleobiodb.org) and [Macrostrat Database](macrostrat.org) tech development initiatives at the University of Wisconsin - Madison.

## Contents
+ [Creative Commons License](#creative-commons-license) # The creative commons license for modules in this repository.
+ [Version and Change Log](#version-and-change-log) # Information about changes to this repository.
+ [communityMatrix.R](#communitymatrixr) # Downloading and cleaning [Paleobiology Database](paleobiodb.org) (PBDB) data, and making a community matrix.
+ [cullMatrix.R](#cullmatrixr) # Culling a communty matrix of depauperate samples and rare taxa.
+ [subsampleRichness.R](#subsamplerichnessr) # A set of subsampling functions for standardizing sampled taxonomic richness.

## Creative Commons License
All code within the [paleobiologyDatabase.R](https://github.com/aazaff/paleobiologyDatabase.R) repository is covered under a Creative Commons License [(CC BY-NC 4.0)](http://creativecommons.org/licenses/by-nc/4.0/). This license requires attribution to Andrew A. Zaffos and Steven M. Holland, and does not allow for commerical use.

## Version and Change Log
This is v0.02 of the paleobiologyDatabase.R repository. The repository has three functional modules: [communityMatrix.R](#communitymatrixr), [cullMatrix.R](#cullmatrixr), and [subsampleRichness.R](#subsamplerichnessr).

+ v.0.02 - Added [subsampleRichness.R](#subsamplerichnessr) module. Changed repository name from CleaningPBDB to paleobiologyDatabase.R. Added new function, softCull(), function to cullMatrix.R module.
+ v.0.01 - Added [communityMatrix.R](#communitymatrixr) and [cullMatrix.R](#cullmatrixr) modules.

## communityMatrix.R
A set of functions written by [Andrew A. Zaffos](https://macrostrat.org/) for downloading data from the Paleobiology Database, and organizing it into a community matrix of samples by taxa, where "samples" are based on  one of the variables in the Paleobiology Database - e.g., Plate ID, Geologic Age.

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

# Parameter DataPBBDB is a dataset downloaded from the PBDB - i.e., using downloadPBDB( )
# Parameter Timescale is a dataset downloaded from Macrostrat - i.e., using downloadTime( )

ConstrainedAges<-constrainAges(DataPBDB=DataPBDB,Timescale=Epochs)
````

##### cleanGenus( )
````
# Cleans the genus field of the PBDB data by removing subgenera and NA's.

# Parameter DataPBBDB is a dataset downloaded from the PBDB - i.e., using downloadPBDB( )

CleanedPBDB<-cleanGenus(DataPBDB)
````

##### communityMatrix( )
````
# Create a community matrix of samples v. species, using elements within one of the PBDB columns
# (e.g., geoplate, early_interval) as the definition of a sample.

# Parameter DataPBBDB is a dataset downloaded from the PBDB - i.e., using downloadPBDB( )
# Parameter SampleDefinition is the column name defining samples

CommunityMatrix<-communityMatrix(DataPBDB,SampleDefinition="geoplate")
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

CulledMatrix<-softMatrix(CommunityMatrix,minOccurrences=5,minDiversity=5)
````

## subsampleRichness.R
Functions for standardizing tasonomic richness.

Can be accessed directly in R using:

````
source("https://raw.githubusercontent.com/aazaff/paleobiologyDatabase.R/master/subsampleRichness.R")
````

##### subsampleEvenness( )
````
# A function that subsamples richness based on evenness. Often referred to as "Shareholder Quorum Subsampling".
# An optimized version of code written John Alroy, by Steven M. Holland.

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
# A multicore version of subsampleEvenness(). Be warned that multicoreEvenness() is not automatically
# faster than subsampleEvenness(), particularly for a low number of trials. It's use is not recommended
# for small abundance vectors or a small numbers of trials.

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
# A multicore version of subsampleIndividuals(). Be warned that multicoreIndividuals() is not automatically
# faster than subsampleIndividuals(), particularly for a low number of trials. It's use is not recommended
# for small abundance vectors or a small numbers of trials.

# Parameter Abundance is a vector of abundances.
# Parameter Quota is the number of individuals to be subsampled. If the Quota is greater than
# the number of individuals, the function will print a warning and return the unstandardized richness.
# Parameter Trials determines how many iterations of the bootstrap are performed, default = 1000
# Parameter Cores sets the number of processor cores, default = 4.

SubsampledRichness<-multicoreIndividuals(Abundance,Quota,Trials=1000,Cores=4)
````
