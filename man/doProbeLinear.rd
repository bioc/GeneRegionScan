
\name{doProbeLinear}
\alias{doProbeLinear}

\title{Do Probe Linear}
\description{
Internal function used to calculate linear correlation between probe level intensity and factorial pData entries such as genotypes.
}
\usage{
    doProbeLinear(object, label, testType="linear model")
}
\arguments{
  \item{object}{A ProbeLevelSet object or a regular ExpressionSet object.}
  \item{label}{An optional character string specifying a column name in the pData of the object. If this argument is given, the gene plot will be colour coded based on the different groups (factors) in the pData entry. If a summaryType other than 'dots' is selected the summarisation is done stratified by the different groups in the pData.}    
  \item{testType}{Character string, defining a statistic procedure to identify especially interesting probes.}
}
\value{A correlation list with p-value and fold-change in columns and probe IDs as rownames. Used by \link{plotStatistics} when decorating plots from \link{plotOnGene}.}
\details{This function is not meant to be called by the user. It does the statistical calculations when called by
\link{plotStatistics}.

One important note about the assumption that the alphabetical order is the same as value order. This usually works with
SNP datatypes of the type AA, AB, BB. However the sorting is locale dependent and can produce odd results, which is why the
value assignments are always printed. See \link[base]{Comparison} for details on this. Further user-control over this aspect of 
the plotting function was omitted, giving more weight to the simplicity of the function interface. In case of problems it 
is recommended to rename the entries in the label so that they will be correctly sorted. 
} 
 
\author{Lasse Folkersen}
\seealso{\code{\link{plotStatistics}}, \code{\link{plotOnGene}}}
\examples{
	data(exampleProbeLevelSet)
	colnames(pData(exampleProbeLevelSet))
	linearModelData<-GeneRegionScan:::doProbeLinear(exampleProbeLevelSet,"genotype3")
	colnames(linearModelData)
	#P value for the probe with ID 679083 in relation to "genotype3"
	print(linearModelData["679083","p-value"])


}

\keyword{documentation}
\keyword{utilities}

