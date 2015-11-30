# CleaningPBDB
Functions for download, cleaning, or culling PBDB data

## communityMatrix
A set of functions written by [Andrew A. Zaffos](https://macrostrat.org/) for downloading data from the Paleobiology Database, and organizing it into a community matrix of samples by taxa, where "samples" are based on  one of the variables in the Paleobiology Database - i.e., Plate ID, Geologic Age.

Can be accessed directly in R using:

````
source("https://raw.githubusercontent.com/aazaff/CleaningPBDB/master/communityMatrix.R")
````

#### downloadPBDB( )
````
# Download data from PBDB by taxonomic group and geologic interval
# Parameter Taxa must be a vector of one or more taxon names (as a character string), no default.
# Parameter StartInterval must be a single interval name accepted by the PBDB, default is "Pliocene"
# Parameter StopInterval must be a single interval name accepted by the PBDB, default is "Pleistocene" 

DataPBDB<-downloadPBDB(Taxa=c("Bivalvia","Gastropoda"),StartInterval="Cambrian",StopInterval="Pleistocene")
````

#### downloadTime( )
````
# Download timescale definitions from Macrostrat
# Parameter Timescale must be a timescale recognized by the macrostrat API, no default
# A list of timescale defs can be seen here https://macrostrat.org/api/defs/timescales?all

Epochs<-downloadTime(TimeScale="international epochs")
````

#### constrainAges( )
````
# Assign fossil occurrences to different ages
# Then remove occurrences that are not temporally constrained to a single interval
# Parameter DataPBBDB is a dataset downloaded from the PBDB - i.e., using downloadPBDB( )
# Parameter Timescale is a dataset downloaded from Macrostrat - i.e., using downloadTime( )

ConstrainedAges<-constrainAges(DataPBDB=DataPBDB,Timescale=Epochs)
````

#### cleanGenus( )
````
# Cleans the genus field of the PBDB data by removing subgenera and NA's
# Parameter DataPBBDB is a dataset downloaded from the PBDB - i.e., using downloadPBDB( )

CleanedPBDB<-cleanGenus(DataPBDB)
````

#### communityMatrix( )
````
# Create a community matrix of samples v. species, using elements within one of the PBDB columns
# (e.g., geoplate, early_interval) as the definition of a sample
# Parameter DataPBBDB is a dataset downloaded from the PBDB - i.e., using downloadPBDB( )
# Parameter SampleDefinition is the column name defining samples

CommunityMatrix<-communityMatrix(DataPBDB,SampleDefinition="geoplate")
````

## cullMatrix
A set of functions written by [Steven M. Holland](http://strata.uga.edu/) for removing depauperate and rare taxa from community matrices of samples by taxa.

Can be accessed directly in R using:

````
source("https://raw.githubusercontent.com/aazaff/CleaningPBDB/master/cullMatrix.R")
````
#### cullMatrix( )
````
# Cull a community matrix of depauperate samples and rare taxa.
# Parameter x is a community matrix, no default.
# Parameter minOccurrences is the minimum number of occurrences for each taxon, default = 2
# Parameter minDiversity is the minimum number of taxa within each sample, default=2

CulledMatrix<-cullMatrix(CommunityMatrix,minOccurrences=5,minDiversity=5)
````

