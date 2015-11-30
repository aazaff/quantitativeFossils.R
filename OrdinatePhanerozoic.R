# Custom functions are camelCase. Arrays, parameters, and arguments are PascalCase
# Dependency functions are not embedded in master functions

# Many functions are repeated between FIGURE modules. This is done so that each figure can be
# generated independently of the other figure scripts - i.e., by only using the code within that module.


######################################### Load Required Libraries ###########################################
# Load Required Libraries
library(plyr)
library(RCurl)
library(doMC)
library(pbapply)
library(vegan)

#############################################################################################################
###################################### FOSSIL DATA Zaffos et al. 2015, ALICE ################################
#############################################################################################################
# Download data from PBDB by taxonomic group and geologic interval=eriod
# Object Taxa must be a vector of taxa
downloadPBDB<-function(Taxa,StartInterval="Pliocene",StopInterval="Pleistocene") {
	Taxa<-paste(Taxa,collapse=",")
	URL<-paste("http://paleobiodb.org/data1.2/occs/list.csv?base_name=",Taxa,"&interval=",StartInterval,",",StopInterval,"&show=paleoloc,phylo&limit=all",sep="")
	File<-read.csv(URL,header=T)
	return(File)
	}

# Download timescales from Macrostrat
downloadTime<-function(Timescale) {
	URL<-paste("https://dev.macrostrat.org/api/defs/intervals?format=csv&timescale=",Timescale,sep="")
	GotURL<-getURL(URL)
	Intervals<-read.csv(text=GotURL,header=T)
	Midpoint<-apply(Intervals[,c("t_age","b_age")],1,median)
	Intervals<-cbind(Intervals,Midpoint)
	return(Intervals)
	}

# Remove subgenera
cleanSubgenera<-function(DataPBDB) {
	SpaceSeparated<-sapply(as.character(DataPBDB[,"genus"]),strsplit," ")
	DataPBDB[,"genus"]<-sapply(SpaceSeparated,function(S) S[1])
	return(DataPBDB)
	}

# Assign fossil occurrences to different ages
sortAges<-function(DataPBDB,TimeScale) {
	DataPBDB[,"early_interval"]<-as.character(DataPBDB[,"early_interval"])
	DataPBDB[,"late_interval"]<-as.character(DataPBDB[,"late_interval"])
	for (i in 1:nrow(TimeScale)) {
		EarlyPos<-which(DataPBDB[,"max_ma"]>TimeScale[i,"t_age"] & DataPBDB[,"max_ma"]<=TimeScale[i,"b_age"])
		DataPBDB[EarlyPos,"early_interval"]<-as.character(TimeScale[i,"name"])
		LatePos<-which(DataPBDB[,"min_ma"]>=TimeScale[i,"t_age"] & DataPBDB[,"min_ma"]<TimeScale[i,"b_age"])
		DataPBDB[LatePos,"late_interval"]<-as.character(TimeScale[i,"name"])
		}
	return(DataPBDB)
	}

# Phanerozoic Community Matrix
phanerozoicMatrix<-function(MarineAges) {
	FinalMatrix<-matrix(0,nrow=length(unique(MarineAges[,"early_interval"])),ncol=length(unique(MarineAges[,"genus"])))
	rownames(FinalMatrix)<-unique(MarineAges[,"early_interval"])
	colnames(FinalMatrix)<-unique(MarineAges[,"genus"])
	for (i in 1:nrow(FinalMatrix)) {
		PlateSubset<-subset(MarineAges,MarineAges[,"early_interval"]==rownames(FinalMatrix)[i])
		ColumnPosition<-match(PlateSubset[,"genus"],colnames(FinalMatrix))
		FinalMatrix[i,ColumnPosition]<-1
		}
	return(FinalMatrix)
	}

# Steve Holland's Culling Functions
cullMatrix <- function(x, minOccurrences=2, minDiversity=2) {
	y <- x
	while (diversityFlag(y, minDiversity) | occurrencesFlag(y, minOccurrences)) {
		y <- cullTaxa(y, minOccurrences)
		y <- cullSamples(y, minDiversity)
	}
	y
}

cullTaxa <- function(x, minOccurrences) {
	PA <- x
 	PA[PA>0] <- 1
 	occurrences <- apply(PA, MARGIN=2, FUN=sum)
 	aboveMinimumOccurrences <- occurrences >= minOccurrences
 	y <- x[ ,aboveMinimumOccurrences]
 	if (length(y)==0) {print("no taxa left!")}
 	y
}

cullSamples <- function(x, minDiversity) {
	PA <- x
 	PA[PA>0] <- 1
 	diversity <- apply(PA, MARGIN=1, FUN=sum)
 	aboveMinimumDiversity <- diversity >= minDiversity
 	y <- x[aboveMinimumDiversity, ]
 	if (length(y[,1])==0) {print("no samples left!")}
 	y
}

occurrencesFlag <- function(x, minOccurrences) {
  	PA <- x
 	PA[PA>0] <- 1
 	occurrences <- apply(PA, MARGIN=2, FUN=sum)
 	if (min(occurrences) < minOccurrences) {
 		flag <- 1
 	}
 	else {
	 	flag <- 0
 	}
 	flag
}

diversityFlag <- function(x, minDiversity) {
 	PA <- x
 	PA[PA>0] <- 1
 	diversity <- apply(PA, MARGIN=1, FUN=sum)
 	if (min(diversity) < minDiversity) {
 		flag <- 1
 	}
 	else {
 		flag <- 0
 	}
 	flag
}

##################################### FOSSILS: Clean, Load, and Prepare Data ################################
# Download from the API
options(timeout=300)
CanonicalTaxa<-c("Bivalvia","Gastropoda","Anthozoa","Brachiopoda","Trilobita","Bryozoa","Nautiloidea","Ammonoidea","Echinoidea")
MarinePhanerozoic<-downloadPBDB(CanonicalTaxa,"Cambrian","Pleistocene")
MarinePhanerozoic<-subset(MarinePhanerozoic,MarinePhanerozoic[,"genus"]!="") # Remove NA Genera
MarinePhanerozoic<-subset(MarinePhanerozoic,is.na(MarinePhanerozoic[,"genus"])!=TRUE) # Remove NA Genera
MarinePhanerozoic<-MarinePhanerozoic[,c("genus","early_interval","late_interval","max_ma","min_ma","paleolat","paleolng","geoplate")]

# Download epoch timescale
Epochs<-downloadTime("international%20epochs")
Epochs<-subset(Epochs,Epochs[,"b_age"]<=541)
rownames(Epochs)<-Epochs[,"name"]

# Sort Ages
MarineAges<-sortAges(MarinePhanerozoic,Epochs)
MarineAges<-MarineAges[MarineAges[,"early_interval"]==MarineAges[,"late_interval"],]
# Remove Terreneuvian because of small sample size and weak GPlates constraint
MarineAges<-subset(MarineAges,MarineAges[,"early_interval"]!="Terreneuvian")

# Clean subgenera
MarineAges<-cleanSubgenera(MarineAges)

# Create PhanerozoicMatrix
PhanerozoicMatrix<-phanerozoicMatrix(MarineAges)

# Cull depauperates
PhanerozoicMatrix<-cullMatrix(PhanerozoicMatrix)

# Perform Ordinations
PhanerozoicDCA<-decorana(PhanerozoicMatrix)
PhanerozoicScores<-scores(PhanerozoicDCA,display="sites")

# Bind output
PhanerozoicScores<-transform(merge(as.data.frame(PhanerozoicScores),Epochs,by="row.names",all=FALSE),row.names=Row.names,Row.names=NULL)

# Plot Result
plot(x=PhanerozoicScores[,1],y=PhanerozoicScores[,2],pch=16,cex=3,col=as.character(PhanerozoicScores[,"color"]),xlab="DCA Axis 1 (\"Time\")",ylab="DCA Axis 2 (\"Unknown\")",las=1)