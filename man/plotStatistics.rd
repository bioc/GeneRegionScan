
\name{plotStatistics}
\alias{plotStatistics}

\title{Plot Statistics}
\description{
Internal function that will assist \link{plotOnGene} in identifying probes that are significantly changed in comparison with a given pData column.
}
\usage{
    plotStatistics(object, probeData, label, summaryType, testType, interval=NULL,
    forcePValue=FALSE, verbose=TRUE, positionVector, ylim, xlim, cutoff=0.2)
}
\arguments{
  \item{object}{A ProbeLevelSet object or a regular ExpressionSet object (in which case a probeData argument is required.}
  \item{probeData}{Optional if a ProbeLevelSet is submitted as object argument. Otherwise it must be a data frame with rownames corresponding to the featureNames of the ExpressionSet and a column named "sequence" with the probe sequences as character strings}
  \item{label}{An optional character string specifying a column name in the pData of the object. If this argument is given, the gene plot will be colour coded based on the different groups (factors) in the pData entry. If a summaryType other than 'dots' is selected the summarisation is done stratified by the different groups in the pData.}    
  \item{summaryType}{Character string specifying one of the following summary methods: 'median', 'mean', 'quartiles' or 'dots' (i.e. no summary). Specifies how all the sample values or all the samples values in a group if 'label' is given, should be summarised. Defaults to 'median'.}
  \item{testType}{Character string, defining a statistic procedure to identify especially interesting probes.}
  \item{interval}{Optional vector of two integers of bp positions. If given, the plot will only include the sequence from gene in the given interval. The x-axis annotation is preserved from original, so this is useful for zooming on specific regions.}
  \item{forcePValue}{Logical. Is used if the testType argument is used. If TRUE all significantly changed probes have P-value given on the plot. If FALSE, only plots with less than 10 probes significant write P-values. Plots can become very cluttered with data if set to TRUE}
  \item{verbose}{TRUE or FALSE}
  \item{positionVector}{A vector of integers and names specifying the position along the x-axis of each probe.}
  \item{ylim}{A vector of two integers. Specifies the decided ylim value passed into the plotting function.}
  \item{xlim}{A vector of two integers. Specifies the decided xlim value passed into the plotting function.}
  \item{cutoff}{Integer specifying at what p-value probes should be circled when using the 'testType' variable. Defaults to 0.2. For cutoffs higher than 0.05, all probes with P >0.05 will be circled in grey instead of black.}
  
    
}
\value{No value, but decorates an existing plotOnGene plot with circles and p-values for all probes that are significant when investigating label with testType}
\details{This function is not meant to be called by the user. It acts as an adapter between \link{plotOnGene} and
either \link{doProbeLinear} or \link{doProbeTTest} depending on the testType variable. 
} 
 
\author{Lasse Folkersen}
\seealso{\code{\link{doProbeLinear}}, \code{\link{doProbeTTest}}, \code{\link{plotOnGene}}}

\keyword{documentation}
\keyword{utilities}

