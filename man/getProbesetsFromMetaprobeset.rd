
\name{getProbesetsFromMetaprobeset}
\alias{getProbesetsFromMetaprobeset}

\title{Get Probeset IDs From metaprobeset IDs}
\description{
Function that return the probesets mapping to a given set of metaprobesets.
}
\usage{
    getProbesetsFromMetaprobeset(annotation, metaprobesets, pythonPath=NULL, mpsToPsFile=NULL)
}
\arguments{
  \item{annotation}{A characther string giving the type of array for which probesets are needed. Will be used to load the .db file from bioconductor. Optional if mps to ps and transcript cluster file is given.}
  \item{metaprobesets}{A vector of characters giving the metaprobeset IDs for which to find the probe set IDs.}
  \item{pythonPath}{Optional character string with the path for Python software. This is only needed for exon arrays. If Python is in path this will be recognised automatically. Python can be downloaded from http://www.python.org.}
  \item{mpsToPsFile}{The location of a transcript cluster file such as HuEx-1\_0-st-v2.r2.dt1.hg18.full.mps. These can be downloaded from http://www.affymetrix.com. Best to use the "full" type file.}
}
\value{A list of all probesets found in the given metaprobesets.}
\details{This function will parse the relevant library
files from the http://www.affymetrix.com website using a python script. Admittedly this is not optional from an R-only
point-of-view, but it works and its fast. The function will be updated when more R-centric ways of parsing exon array 
annotation data is available. 

Alternatively this data can just as well be retrieved from the web, but in some cases this function is faster an easier. 
}
\author{Lasse Folkersen}
\seealso{\code{\link{getLocalProbeIntensities}}, \code{\link{getProbesetsFromRegionOfInterest}}, \code{\link{getMetaprobesetsFromRegionOfInterest}}}
\examples{

	\dontrun{
	#must supply mpsToPsFile and transcriptClustersFile for this to work
	probesets<-getProbesetsFromMetaprobeset("notusedhere", c("3218528","2423669"), transcriptClustersFile=transcriptClustersFile, mpsToPsFile=mpsToPsFile)
	}
}
\keyword{documentation}
\keyword{utilities}

