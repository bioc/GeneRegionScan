
\name{getLocalProbeIntensities}
\alias{getLocalProbeIntensities}

\title{Get Probe Intensities locally}
\description{
Function that will create a ProbeLevelSet from cel files.
}
\usage{
    getLocalProbeIntensities(listOfProbesets, celfilePath, annotation=NULL, aptCelExtractPath=NULL, pgfPath=NULL, clfPath=NULL,
    cdfPath=NULL, verbose=TRUE)
}
\arguments{
  \item{listOfProbesets}{Either a character vector with the names of the probesets from which probes should be included in the ProbeLevelSet, or the path name of a file containing this}
  \item{celfilePath}{The path to a folder that contains the cel files of interest.}
  \item{annotation}{Character string specifying which type of arrays is investigated.}
  \item{aptCelExtractPath}{The path to the apt-cel-extract file from the Affymetrix Power Tools package (including the filename itself). You can try to leave it as NULL and see what happens. If the tool is in path, it will work anyway. If you have specified the argument before it will be remembered. If you are using an OS for which executables have been included in the package (win32 and linux64), it will use that.}
  \item{pgfPath}{The path to a pgf file for the exon array of interest. This argument is mutually exclusive with the cdfPath argument. These files can be downloaded from www.affymetrix.com for the array of interest. Once given for a particular annotation, the location is saved for future use and can be given as NULL next time.}
  \item{clfPath}{The path to a clf file for the exon array of interest. This argument is mutually exclusive with the cdfPath argument. These files can be downloaded from www.affymetrix.com for the array of interest. Once given for a particular annotation, the location is saved for future use and can be given as NULL next time.}
  \item{cdfPath}{The path to a cdf file for the array of interest. This argument is mutually exclusive with the pgfPath and clfPath arguments. These files can be downloaded from www.affymetrix.com for the array of interest.  Once given for a particular annotation, the location is saved for future use and can be given as NULL next time.}
  \item{verbose}{TRUE or FALSE.}
}
\value{A ProbeLevelSet with probe intensity values for all probes in the specified probesets. The ProbeLevelSet 
inherits the ExpressionSet, but in addition it has the sequence of each probe stored in the featureData of the set.}
\details{This is the workhorse function for generating ProbeLevelSets. It has two main functions:
1) It extracts the probe level intensity, normalizes it to all probes in the data using quantiles normalization and transfers 
the results to the exprs part of a ProbeLevelSet. The extraction and normalization parts are done by a call to 
the Affymetrix Power Tools (APT) function apt-cel-extract. The analysis switch for this call is '-a quant-norm,pm-only'
and the documentation for APT can be used for further reference. 
2) It extracts the probe sequences of all probes of interest. This is done using the [arrayname]probe packages for the 3'IVT type arrays
for which they are available. For the exon type arrays, for which this is not available it is done by direct parsing of the pgf
files using the readPgf function of the affxparser package. The results of this is saved in the featureData of the ProbeLevelSet}
\author{Lasse Folkersen}
\seealso{\code{\link{getServerProbeIntensities}}, \code{\link{plotOnGene}}}
\examples{

	\dontrun{
	#must correct paths and give cel files before this example will work
	listOfProbesets<-c("10321_at")
	celfilePath<-path_to_some_cel_files
	aptCelExtractPath<-"~/apt/apt-cel-extract"
	cdfPath<-"~/hgu133plus2.cdf"
	getLocalProbeIntensities(listOfProbesets, celfilePath, annotation="hgu133plus2", aptCelExtractPath=aptCelExtractPath, cdfPath=cdfPath)
	}
}
\keyword{documentation}
\keyword{utilities}


