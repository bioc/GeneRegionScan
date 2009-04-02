
\name{getServerProbeIntensities}
\alias{getServerProbeIntensities}

\title{Get Probe Intensities from a server}
\description{
Wrapper around \link{getLocalProbeIntensities} that is designed to be run easily on remote computers.
}
\usage{
    getServerProbeIntensities(listOfProbesets, celfilePath, annotation=NULL, aptCelExtractPath=NULL, pgfPath=NULL, clfPath=NULL, cdfPath=NULL, serveradress="localhost", username=NULL, password=NULL, plinkPath=NULL, pscpPath=NULL, verbose=TRUE)}
\arguments{
  \item{listOfProbesets}{Either a character vector with the names of the probesets from which probes should be included in the ProbeLevelSet, or the path name of a file containing this}
  \item{celfilePath}{The path to a folder that contains the cel files of interest. If serveradress is different from "localhost", this should be the path on the remote computer}
  \item{annotation}{Character string specifying which type of arrays is investigated. If a pgf file and clf file is specified this is optional.}
  \item{aptCelExtractPath}{The path to the apt-cel-extract file from the Affymetrix Power Tools package (including the filename itself). You can try to leave it as NULL and see what happens. If the tool is in path, it will work anyway. If you have specified the argument before it will be remembered. If you are using an OS for which executables have been included in the package (win32 and linux64), it will default to using that.}
  \item{pgfPath}{The path to a pgf file for the exon array of interest. This argument is mutually exclusive with the cdfPath argument. These files can be downloaded from www.affymetrix.com for the array of interest. Once given for a particular annotation, the location is saved for future use and can be given as NULL next time.}
  \item{clfPath}{The path to a clf file for the exon array of interest. This argument is mutually exclusive with the cdfPath argument. These files can be downloaded from www.affymetrix.com for the array of interest. Once given for a particular annotation, the location is saved for future use and can be given as NULL next time.}
  \item{cdfPath}{The path to a cdf file for the array of interest. This argument is mutually exclusive with the pgfPath and clfPath arguments. These files can be downloaded from www.affymetrix.com for the array of interest. Once given for a particular annotation, the location is saved for future use and can be given as NULL next time.}
  \item{serveradress}{Character string with the IP adress to a remote server. Defaults to "localhost", in which case the algorithm is run locally.}
  \item{username}{The username for the remote server. Optional if serveradress is "localhost".}
  \item{password}{The password for the remote server. Optional if serveradress is "localhost".}
  \item{plinkPath}{The path to the plink program. Can be found at http://www.chiark.greenend.org.uk/~sgtatham/putty/download.html. Used to transfer commands to the remote server. Optional if serveradress is "localhost".}
  \item{pscpPath}{The path to the pscp program. Can be found at http://www.chiark.greenend.org.uk/~sgtatham/putty/download.html. Used to transfer commands to the remote server. Optional if serveradress is "localhost".}
  \item{verbose}{TRUE or FALSE.}
}
\value{A ProbeLevelSet with probe quantiles normalized intensity values for all probes in the specified probesets. The ProbeLevelSet inherits the ExpressionSet, but in addition it has the sequence of each probe stored in the featureData of the set.}
\details{This function is a wrapper around \link{getLocalProbeIntensities}. The main function of the getServerProbeIntensities 
function is to save all the data needed to start a getLocalProbeIntensities run, send it to a remote server, start the run, and return the data.
This is useful if you have access to a fast server and want to work with large datasets of exon arrays.}
\author{Lasse Folkersen}
\seealso{\code{\link{getLocalProbeIntensities}}, \code{\link{plotOnGene}}}
\examples{

	\dontrun{
	#must correct paths, username and password before this example will work
	listOfProbesets<-c("10321_at")
	celfilePath<-path_to_some_cel_files #remote
	aptCelExtractPath<-"~/apt/apt-cel-extract" #remote
	cdfPath<-"~/hgu133plus2.cdf" #remote
	username<-"user1"
	password<-"ZaphodBeeblebrox"
	plinkPath<-"~/plink" #local
	pscpPath<-"~/pscp" #local
	getServerProbeIntensities(listOfProbesets, celfilePath, annotation="hgu133plus2", aptCelExtractPath=aptCelExtractPath, cdfPath=cdfPath,
	username=username, password=password, plinkPath=plinkPath, pscpPath=pscpPath)

	}
}
\keyword{documentation}
\keyword{utilities}


