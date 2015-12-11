# Custom functions are camelCase. Arrays, parameters, and arguments are PascalCase
# Dependency functions are not embedded in master functions

######################################### Load Required Libraries ###########################################
# Load Required Libraries
if (require("RCurl")==FALSE) {
	install.packages("RCurl")
	library("RCurl")
	} else {
		library("RCurl")
		}

downloadPBDB<-function(Taxa,StartInterval="Pliocene",StopInterval="Pleistocene") {
	Taxa<-paste(Taxa,collapse=",")
	URL<-paste("https://paleobiodb.org/data1.2/occs/list.csv?base_name=",Taxa,"&interval=",StartInterval,",",StopInterval,"&show=paleoloc,phylo&limit=all",sep="")
	GotURL<-getURL(URL)
	File<-read.csv(text=GotURL,header=T)
	return(File)
	}

# Download timescales from Macrostrat
downloadTime<-function(Timescale) {
	Timescale<-gsub(" ","%20",Timescale)
	URL<-paste("https://dev.macrostrat.org/api/defs/intervals?format=csv&timescale=",Timescale,sep="")
	GotURL<-getURL(URL)
	Intervals<-read.csv(text=GotURL,header=T)
	Midpoint<-apply(Intervals[,c("t_age","b_age")],1,median)
	Intervals<-cbind(Intervals,Midpoint)
	return(Intervals)
	}

# Remove subgenera and NAs
cleanGenus<-function(DataPBDB) {
	DataPBDB<-subset(DataPBDB,DataPBDB[,"genus"]!="") # Remove NA Genera
	DataPBDB<-subset(DataPBDB,is.na(DataPBDB[,"genus"])!=TRUE) # Remove NA Genera
	SpaceSeparated<-sapply(as.character(DataPBDB[,"genus"]),strsplit," ")
	DataPBDB[,"genus"]<-sapply(SpaceSeparated,function(S) S[1])
	return(DataPBDB)
	}

# Assign fossil occurrences to different ages
# Then remove occurrences that are not temporally constrained to a single interval
constrainAges<-function(DataPBDB,Timescale) {
	DataPBDB[,"early_interval"]<-as.character(DataPBDB[,"early_interval"])
	DataPBDB[,"late_interval"]<-as.character(DataPBDB[,"late_interval"])
	for (i in 1:nrow(Timescale)) {
		EarlyPos<-which(DataPBDB[,"max_ma"]>Timescale[i,"t_age"] & DataPBDB[,"max_ma"]<=Timescale[i,"b_age"])
		DataPBDB[EarlyPos,"early_interval"]<-as.character(Timescale[i,"name"])
		LatePos<-which(DataPBDB[,"min_ma"]>=Timescale[i,"t_age"] & DataPBDB[,"min_ma"]<Timescale[i,"b_age"])
		DataPBDB[LatePos,"late_interval"]<-as.character(Timescale[i,"name"])
		}
	DataPBDB<-DataPBDB[DataPBDB[,"early_interval"]==DataPBDB[,"late_interval"],] # Remove taxa that range through
	return(DataPBDB)
	}

# Create a community matrix of samples v. species, using elements within one of the PBDB columns
# (e.g., geoplate, early_interval) as the definition of a sample. This is a presence-absence matrix.
presenceMatrix<-function(DataPBDB,SampleDefinition="geoplate") {
	FinalMatrix<-matrix(0,nrow=length(unique(DataPBDB[,SampleDefinition])),ncol=length(unique(DataPBDB[,"genus"])))
	rownames(FinalMatrix)<-unique(DataPBDB[,SampleDefinition])
	colnames(FinalMatrix)<-unique(DataPBDB[,"genus"])
	for (i in 1:nrow(FinalMatrix)) {
		SampleSubset<-subset(DataPBDB,DataPBDB[,SampleDefinition]==rownames(FinalMatrix)[i])
		ColumnPosition<-match(SampleSubset[,"genus"],colnames(FinalMatrix))
		FinalMatrix[i,ColumnPosition]<-1
		}
	return(FinalMatrix)
	}

# Create a community matrix of samples v. species, using elements within one of the PBDB columns
# (e.g., geoplate, early_interval) as the definition of a sample. This is an "abundance" matrix which uses
# the number of occurrences.
abundanceMatrix<-function(DataPBDB,SampleDefinition="geoplate") {
	FinalMatrix<-matrix(0,nrow=length(unique(DataPBDB[,SampleDefinition])),ncol=length(unique(DataPBDB[,"genus"])))
	rownames(FinalMatrix)<-unique(DataPBDB[,SampleDefinition])
	colnames(FinalMatrix)<-unique(DataPBDB[,"genus"])
	for (i in 1:nrow(FinalMatrix)) {
		SampleSubset<-subset(DataPBDB,DataPBDB[,SampleDefinition]==rownames(FinalMatrix)[i])
		Abundances<-table(SampleSubset[,"genus"])
		Abundances<-Abundances[Abundances>0]
		ColumnPosition<-match(names(Abundances),colnames(FinalMatrix))
		FinalMatrix[i,ColumnPosition]<-Abundances
		}
	return(FinalMatrix)
	}
