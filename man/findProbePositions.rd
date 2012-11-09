
\name{findProbePositions}
\alias{findProbePositions}
\alias{findProbePositions,ExpressionSet-method}

\title{Find positions of probes on a gene}
\description{
Function that will the location of probes in a gene based on their sequence.
}
\usage{
    findProbePositions(object, gene, probeData=NULL,interval=NULL,directions="all",verbose=TRUE)
}
\arguments{
  \item{object}{A ProbeLevelSet object or a regular ExpressionSet object (in which case a probeData argument is required). See \link{getLocalProbeIntensities} and related functions on how to create a ProbeLevelSet.}
  \item{gene}{A number of gene sequences as DNAstring, DNAStringsSets, or character-vectors with sequence.}
  \item{probeData}{Optional if a ProbeLevelSet is submitted as object argument. Otherwise, it must be a data frame with rownames corresponding to the featureNames of the ExpressionSet and a column named "sequence" with the probe sequences as character strings}
  \item{interval}{Optional vector of two integers of bp positions. If given, the plot will only include the sequence from gene in the given interval. The x-axis annotation is preserved from original, so this is useful for zooming on specific regions.}
  \item{verbose}{TRUE or FALSE}
  \item{directions}{A character vector of the matching-directions that should be scanned (which combinations of complementary and reverse). Defaults to "all" which is shorthand for all possible directions, but can take anything from: c("matchForwardSense", "matchForwardAntisense", "matchReverseSense", "matchReverseAntisense")}
  
    
}
\value{A vector of positions of each probe in the object ProbeLevelSet, with names being probe ids}
\details{
	This function is principally used by the \code{\link{plotOnGene}} to assign positions of each probe relative to the gene sequence
	of interest. In a recent version of GeneRegionScan it was separated as a discrete function because of its use in alternative plotting
}
\author{Lasse Folkersen}
\seealso{\code{\link{geneRegionScan}}, \code{\link{plotOnGene}}, \code{\link{plotCoexpression}}}
\examples{

	data(exampleProbeLevelSet)
	
	findProbePositions(exampleProbeLevelSet, mrna)
	
}
\keyword{documentation}
\keyword{utilities}

