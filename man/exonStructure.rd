
\name{exonStructure}
\alias{exonStructure}

\title{Add exon-structure to plots}
\description{
Function that will paint the exon structure of a gene on the plots obtained by plotOnGene.
}
\usage{
    exonStructure(mrna, genome, maxMismatch=4, y=0)
}
\arguments{
  \item{mrna}{A gene sequence formatted as DNAstring, character-vectors or \link[Biostrings]{readFASTA} output.}
  \item{genome}{A number of gene sequences as DNAstring, vectors of DNAStrings, character-vectors or \link[Biostrings]{readFASTA} outputs.}
  \item{maxMismatch}{Integer. The maximum number of mismatches per exon that can be allowed before the exon is not allocated at the position in the template mrna. Defaults to 4.}
  \item{y}{Numeric. The vertical position of the exon structure (if ylim is changed)}    
}
\value{No value, but plots the layout of exons in a gene on the product of a call to \link{plotOnGene}.}
\details{When given a sequence of the DNA divided by exons and a sequence of the corresponding mRNA string, this function
will plot the layout of exons along the length of the x-axis on the current device. The sequences have to be given
as FASTA sequences produced by the \link[Biostrings]{readFASTA} function in the Biostrings package. Furthermore the genome
must be divided with an entry for each exon. This is easily done by downloading the genome sequences of the gene-of-interest
from http://genome.ucsc.edu and specifying "One FASTA record per region".

The function works by sequentially comparing each exon to the mrna. The location of the match start and end is taken
as the exon boundaries and plotted. If more than one match is found a warning is given. If no match is found for an exon
this is printed, but otherwise ignored. The mrna can in fact also be DNA sequence with introns. The important thing 
is that it serves as template for the exon matching.

Importantly, the exon-numbers technically refer to "number of exon in investigated transcript". For example if 
an DNA with exon structure for an isoform which does not include all exons in the gene is investigated, then there will 
be skips of exon numbers. To avoid this, the DNA for an isoform which do include all exons could be used. 
However, it is really more a biological issue: different sources can differ on where to start 
the counting for a given gene in any case.  
	  
} 
 
\author{Lasse Folkersen}
\seealso{\code{\link{geneRegionScan}}, \code{\link{plotOnGene}}}
\examples{

	data(exampleProbeLevelSet)
	plotOnGene(exampleProbeLevelSet, mrna, label="gender", testType="wilcoxon")
	exonStructure(mrna, genomic, maxMismatch=2)
}
\keyword{documentation}
\keyword{utilities}

