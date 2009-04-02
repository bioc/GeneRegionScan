
\name{translateSampleNames}
\alias{translateSampleNames}
\alias{translateSampleNames,ExpressionSet-method}


\title{Translate sample names using translation file}
\description{
Function that will change sample names in an ExpressionSet using a translation file. 
}
\usage{
    translateSampleNames(object, translationFile, from, to)
}
\arguments{
  \item{object}{An ExpressionSet}
  \item{translationFile}{Character string with the path to a tab-separated text file with translations of names. Alternatively a data frame containing translation information.}
  \item{from}{Character string with the translationFile column name containing the values to translate from.}
  \item{to}{Character string with the translationFile column name containing the values to translate to.}
}
\value{The same ExpressionSet given as argument, but with sampleNames changed as specified in translationFile}
\details{
Function that can translate the sampleNames of an ExpressionSet when given translation information.
The translationFile argument gives the translation between sampleNames as they are and the sampleNames as they should be. 
It can either be a character string with the path to a tab-separated text-file with headers or a directly a data frame with
the necessary information. The variables 'from' and 'to' specify the column names that contains the translation information.
They must be present in the translationFile and the entries in 'from' must be present in the sampleNames of the 
current ExpressionSet. Extra entries in the translationFile are omitted but samples in the ExpressionSet without 
samples in the translationFile will raise an error.

This package is particularly useful after running cel-files through the Affymetrix Power Tools 
(with \link{getLocalProbeIntensities} for example), since this program has a tendency to change any odd 
character in the filename to something else.

}
\author{Lasse Folkersen}

\examples{

	data(exampleProbeLevelSet)
	
	fromThese<-sampleNames(exampleProbeLevelSet)
	
	toThese<-sub(".cel","", sampleNames(exampleProbeLevelSet))
	toThese<-sub("X","", toThese)
	

	translationFile<-cbind(fromThese, toThese)
	translationFile<-as.data.frame(translationFile)
	
	translateSampleNames(exampleProbeLevelSet, translationFile, from="fromThese", to="toThese")
				
	
}
\keyword{documentation}
\keyword{utilities}


