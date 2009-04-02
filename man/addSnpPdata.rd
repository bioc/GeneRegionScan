
\name{addSnpPdata}
\alias{addSnpPdata}
\alias{addSnpPdata,ExpressionSet-method}


\title{Add SNP data from hapmap-style data file to pData of ExpressionSet}
\description{
Function that will add genotypes to an ExpressionSet when given a data file with genotypes in the same format as outputted by http://www.hapmap.org.
}
\usage{
    addSnpPdata(object, listOfSnps)
}
\arguments{
  \item{object}{A ProbeLevelSet object or a regular ExpressionSet object (in which case a probeData argument is required). See \link{getLocalProbeIntensities} and related functions on how to create a ProbeLevelSet.}
  \item{listOfSnps}{A character string giving the path of a file containing SNP data in the same format as used by the export function of www.hapmap.org.}
}
\value{The same ProbeLevelSet or ExpressionSet given as argument, but with the SNP type data added to the pData and with the extra information of the listOfSnps file added to the notes section.}
\details{
A function that takes an ExpressionSet and the path for a SNP-file. The SNP-file should have the same format as 
files downloaded from http://www.hapmap.org. The function then checks if the sampleNames are present in the SNP-file. If they are, it decorates the ExpressionSet with the SNPs.
If the snp name already exists as a column name in the pData, the function tries to merge the new information from listOfSnps with the new. Original NA values are overwritten by new values. Original values that are the same
as new values are kept the same. Original values that are not the same as new values are preserved in their original state and a warning is given.
 
}
\author{Lasse Folkersen}
\seealso{\code{\link{getLocalProbeIntensities}},\code{\link{plotOnGene}}}
\examples{

	\dontrun{
	hapmapformatdata<-"~/somefile.txt"
	probelevelsetwithsnps<-addSnpPdata(probelevelset,hapmapformatdata)
	}
}
\keyword{documentation}
\keyword{utilities}

