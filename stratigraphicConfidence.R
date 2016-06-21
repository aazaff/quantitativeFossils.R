# From Marshall's (1990) adaptation of Strauss and Sadler.
estimateExtinction <- function(OccurrenceAges, ConfidenceLevel=.95)  {
  # Find the number of unique "Horizons"
  NumOccurrences<-length(unique(OccurrenceAges))-1
  Alpha<-((1-ConfidenceLevel)^(-1/NumOccurrences))-1
  Lower<-min(OccurrenceAges)
  Upper<-min(OccurrenceAges)-(Alpha*10)
  return(setNames(c(Lower,Upper),c("Earliest","Latest")))
  }

