
\name{getLocalMetaprobeIntensities}
\alias{getLocalMetaprobeIntensities}

\title{Get Metaprobe Intensities locally}
\description{
Function that will create an expressionset from cel files.
}
\usage{
	getLocalMetaprobeIntensities(celfilePath, analysis="rma", metaProbeSetsFile=NULL, annotation=NULL, aptProbesetSummarizePath=NULL,
	pgfPath=NULL, clfPath=NULL, cdfPath=NULL, verbose=TRUE)

}


\arguments{
  \item{celfilePath}{The path to a folder that contains the cel files of interest.}
  \item{analysis}{The analysis string passed to Affymetrix Power Tools. "RMA" and "RMA-sketch" are recommended starting points. See APT documentation for complete list of possibilities.}
  \item{metaProbeSetsFile}{Path to a file containing meta probe set information. These can be downloaded from www.affymetrix.com and has names like HuEx\-1\_0\-st\-v2.r2.dt1.hg18.full.mps. This choice also decides which subset of the meta probe sets that will be used: Core, Extended or Full.}
  \item{annotation}{Character string specifying which type of arrays is investigated.}
  \item{aptProbesetSummarizePath}{The path to the apt-probeset-summarize file from the Affymetrix Power Tools package (including the filename itself). You can try to leave it as NULL and see what happens. If you have specified the argument before it will be remembered.}
  \item{pgfPath}{The path to a pgf file for the exon array of interest. This argument is mutually exclusive with the cdfPath argument. These files can be downloaded from www.affymetrix.com for the array of interest. Once given for a particular annotation, the location is saved for future use and can be given as NULL next time.}
  \item{clfPath}{The path to a clf file for the exon array of interest. This argument is mutually exclusive with the cdfPath argument. These files can be downloaded from www.affymetrix.com for the array of interest. Once given for a particular annotation, the location is saved for future use and can be given as NULL next time.}
  \item{cdfPath}{The path to a cdf file for the array of interest. This argument is mutually exclusive with the pgfPath and clfPath arguments. These files can be downloaded from www.affymetrix.com for the array of interest.  Once given for a particular annotation, the location is saved for future use and can be given as NULL next time.}
  \item{verbose}{TRUE or FALSE.}
}
\value{An Expressionset with meta probe set values for all meta probes in the specified metaProbeSetFile.}
\details{This function is a simple wrapper around Affymetrix Power Tools (APT). It is useful for importing
meta probe set data from cel files to R environments. It will remember the location of annotation files,
after the first use, which removes some of the typing otherwise included in using APT.
This function has now been updated so that a windows edition of apt-probeset-summarize.exe from APT-1.12.0 is 
included in the installation package, so that pre-processing of cel files can be done with less setup. Linux users
will still have to install APT for this function. 
}
\author{Lasse Folkersen}
\seealso{\code{\link{getLocalProbeIntensities}}, \code{\link{plotOnGene}}}
\examples{

	\dontrun{
	#must correct paths and give cel files before this example will work
	listOfProbesets<-c("10321_at")
	celfilePath<-path_to_some_cel_files
	aptProbesetSummarizePath<-"~/apt/apt-cel-extract"
	cdfPath<-"~/hgu133plus2.cdf"
	getLocalMetaprobeIntensities(celfilePath, annotation="hgu133plus2", aptProbesetSummarizePath=aptProbesetSummarizePath, cdfPath=cdfPath)
	}
}
\keyword{documentation}
\keyword{utilities}


