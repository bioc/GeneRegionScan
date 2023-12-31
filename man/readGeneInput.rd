
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
  \item{gene}{A number of gene sequences as DNAString, DNAStringSets or, character-vectors with sequence.} 
  \item{genename}{Optional character string specifying a gene name to include
    in the plot. If not included and a FASTA sequence is given, it will default
    to the name in the FASTA sequence. Otherwise it will default to 'Unknown
    genename'.} 
  \item{verbose}{TRUE or FALSE.}
}
\value{A list of all sequences and their names, in exactly the same format as obtained with the deprecated function readFASTA.}
\details{This function is not meant to be run directly by the user. It will take a number of genes either as a vector 
of characters, a path to a FASTA format file, as a DNAstrings or DNAStringSet. It will then output 
them in FASTA format for use with other functions. Optional variable genename forces a new name.

The primary objective of this function is to make the sequence input simpler for other functions.
 
}
\author{Lasse Folkersen}
\seealso{\code{\link{geneRegionScan}}, \code{\link{plotOnGene}}}
\examples{

	somegene<-"ATACCTTGTAGGACCTGATGATAGATGCATAGTAATATCGTA"
	genename<-"My favourite gene"
	GeneRegionScan:::readGeneInput(somegene,genename=genename)

}
\keyword{documentation}
\keyword{utilities}




