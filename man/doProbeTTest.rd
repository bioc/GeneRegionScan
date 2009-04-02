
\name{doProbeTTest}
\alias{doProbeTTest}

\title{Do Probe T-Test}
\description{
Internal function used to calculate t-tests between probe level intensity and pData entries with two levels.
}
\usage{
    doProbeTTest(object, label, testType="students")
}
\arguments{
  \item{object}{A ProbeLevelSet object or a regular ExpressionSet object.}
  \item{label}{An optional character string specifying a column name in the pData of the object. If this argument is given, the gene plot will be colour coded based on the different groups (factors) in the pData entry. If a summaryType other than 'dots' is selected the summarisation is done stratified by the different groups in the pData.}    
  \item{testType}{Character string, defining a statistic procedure to identify especially interesting probes.}
}
\value{A correlation list with p-value and fold-change in columns and probe IDs as rownames. Used by \link{plotStatistics} when decorating plots from \link{plotOnGene}.}
\details{This function is not meant to be called by the user. It does the statistical calculations when called by
\link{plotStatistics} 
} 
 
\author{Lasse Folkersen}
\seealso{\code{\link{plotStatistics}}, \code{\link{plotOnGene}}}
\examples{

	data(exampleProbeLevelSet)
	colnames(pData(exampleProbeLevelSet))
	tTestData<-GeneRegionScan:::doProbeTTest(exampleProbeLevelSet,"gender")
	colnames(tTestData)
	
	#P value for the probe with ID 679083 in relation to "gender"
	print(tTestData["679083","p-value"])

}


\keyword{documentation}
\keyword{utilities}

