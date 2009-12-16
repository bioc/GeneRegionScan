
\name{readGeneInput}
\alias{readGeneInput}

\title{Standardize reading of gene inputs}
\description{
Internal function that will standardise the input of the many different form of sequences that can be used in Bioconductor.
}
\usage{
    readGeneInput(gene, genename=NULL, verbose=TRUE)
}
\arguments{
  \item{gene}{}
  \item{genename}{}
  \item{verbose}{}
}
\value{A list of all sequences and their names, in exactly the same format as obtained with the \link[Biostrings]{readFASTA} function of the Biostrings package.}
\details{This function is not meant to be run directly by the user. It will take a number of genes either as a vector 
of characters, a path to a FASTA format file, as a vector of DNAstrings or as a readFASTA format. It will then output 
them in FASTA format for use with other functions. Optional variable genename forces a new name.

The primary objective of this function is to make the sequence input simpler for other functions.
 
}
\author{Lasse Folkersen}
\seealso{\code{\link{geneRegionScan}},\code{\link{plotOnGene}}}
\examples{

	somegene<-"ATACCTTGTAGGACCTGATGATAGATGCATAGTAATATCGTA"
	genename<-"My favourite gene"
	GeneRegionScan:::readGeneInput(somegene,genename=genename)

}
\keyword{documentation}
\keyword{utilities}




