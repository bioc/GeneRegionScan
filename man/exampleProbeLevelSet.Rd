\name{exampleProbeLevelSet}
\docType{data}
\alias{exampleProbeLevelSet}
\title{Example Dataset of class 'ProbeLevelSet'}
\description{
This ProbeLevelSet is derived from Human exon ST 1.0 cel files downloaded from the http://www.ncbi.nlm.nih.gov/geo/.
The pData in the ProbeLevelSet is entirely fictional, and the intensity have been changed for some probes to show
a point about correlation to genotype.
  }
\usage{data(exampleProbeLevelSet)}
\format{
	A \code{\link[GeneRegionScan:ProbeLevelSet-class]{ProbeLevelSet}} object containing probe level intensity data
	for a region on chromosome 2.
}
\examples{
	data(exampleProbeLevelSet)
	pData(exampleProbeLevelSet)
	exprs(exampleProbeLevelSet)[1:10,]
	
}


\author{Lasse Folkersen}
\keyword{datasets}