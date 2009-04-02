
\name{plotCoexpression}
\alias{plotCoexpression}
\alias{plotCoexpression,ExpressionSet-method}


\title{Plot Coexpression of probes in a ProbeLevelSet}
\description{
Function that will investigate all possible pairings in a set of probes, calculate the Pearson correlation coefficient
and plot them in a meaningful way
}
\usage{
    plotCoexpression(object, gene, probeData=NULL, verbose=TRUE, directions="all", correlationCutoff=0.5,
    probeLevelInfo=c("probeid"))
}
\arguments{
  \item{object}{A ProbeLevelSet object or a regular ExpressionSet object (in which case a probeData argument is required). See \link{getLocalProbeIntensities} and related functions on how to create a ProbeLevelSet.}
  \item{gene}{A number of gene sequences as DNAstring, vectors of DNAStrings, character-vectors or \link[Biostrings]{readFASTA} outputs.}
  \item{probeData}{Optional if a ProbeLevelSet is submitted as object argument. Otherwise it must be a data frame with rownames corresponding to the featureNames of the ExpressionSet and a column named "sequence" with the probe sequences as character strings}
  \item{verbose}{TRUE or FALSE}
  \item{directions}{A character vector of the matching-directions that should be scanned (which combinations of complementary and reverse). Defaults to "all" which is shorthand for all possible directions, but can take anything from: c("matchForwardSense", "matchForwardAntisense", "matchReverseSense", "matchReverseAntisense")}
  \item{correlationCutoff}{A number between 0 and 1. The limit at which Pearson correlation (in absolute values) should not be plotted below. Defaults to 0.5}
  \item{probeLevelInfo}{The information about each probe to include in the plot. Should be a vector of one or more of the following elements: probeid, probesetid, sequence. Default is only probeid.}
  
    
}
\value{No value, but plots a hapmap style plot of correlation values between all probes}
\details{This function takes a ProbeLevelSet or an ExpressionSet + probeData and the sequence of a gene. It then calculates
pairwise Pearson correlation coefficients between all possible combinations of probes. Then it assigns each probe to a location 
along the length of the gene and plots a relational graph showing which probes has high correlation coefficients. The correlation coefficients are sorted by
absolute values meaning that it will also include the negative correlations.}
\author{Lasse Folkersen}
\seealso{\code{\link{geneRegionScan}}, \code{\link{plotOnGene}}}
\examples{
	
	data(exampleProbeLevelSet)
    plotCoexpression(exampleProbeLevelSet, mrna, correlationCutoff=0.7, probeLevelInfo=c("probeid","sequence"))

}
\keyword{documentation}
\keyword{utilities}

