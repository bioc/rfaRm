\name{rfamSequenceSearch}
\alias{rfamSequenceSearch}
\title{Performs a sequence search of the Rfam database
}
\description{
Performs a search of the Rfam database by a provided RNA sequence, and retrieves
high-scoring hits of Rfam families with different regions of the provided
sequence.
}
\usage{
rfamSequenceSearch(sequence)
}
\arguments{
  \item{sequence}{string with an RNA sequence to be searched against the Rfam
  database. Should contain only standard RNA symbols (i.e., "A", "U", "G" and
  "C"), and no spaces or newlines.
  }
}
\value{
A nested list where each element of the top-level list represents a high-scoring
hit with the Rfam families. Each of the top-level list elements is a list in
itself, containing the following elements that describe the hit:
  \item{rfamAccession }{Rfam accession of the Rfam family with which a hit was
  found}
  \item{bitScore }{Bit score for the hit with an Rfam family
  }
  \item{eValue }{Expectation value for the with an Rfam family
  }
  \item{alignmentStartPositionQuerySequence }{Start position in the query
  sequence of the sequence region that resulted in the hit
  }
  \item{alignmentStartPositionHitSequence }{Start position in the Rfam family
  consensus sequence region with which the hit was found
  }
  \item{alignmentEndPositionHitSequence }{End position in the Rfam family
  consensus sequence region with which the hit was found}
  \item{alignmentQuerySequence }{Sequence region of the query RNA sequence
  provided to search the Rfam database, aligned with the corresponding region
  of the consensus sequence of the Rfam family with which the hit was found}
  \item{alignmentMatch }{String describing the matches between the aligned
  regions of the query sequence and the consensus sequence of the Rfam family}
  \item{alignmentHitSequence }{Sequence region of the consensus sequence of the
  Rfam family with which the hit was found, aligned to the corresponding region
  of the query RNA sequence}
  \item{alignmentSecondaryStructure }{Secondary structure of the region
  of the consensus sequence of the Rfam family with which the hit was found
  }
}
\references{
Ioanna Kalvari, Joanna Argasinska, Natalia Quinones-Olvera, Eric P Nawrocki,
Elena Rivas, Sean R Eddy, Alex Bateman, Robert D Finn, Anton I Petrov,
Rfam 13.0: shifting to a genome-centric resource for non-coding RNA families,
Nucleic Acids Research, Volume 46, Issue D1, 4 January 2018, Pages D335–D342,
https://doi.org/10.1093/nar/gkx1038

https://docs.rfam.org/en/latest/api.html

https://www.tbi.univie.ac.at/RNA/ViennaRNA/doc/html/rna_structure_notations.html
}
\examples{
# Search the Rfam database for hits with a specific sequence, and store the
# results in a nested list

searchHits <- rfamSequenceSearch("GGAUCUUCGGGGCAGGGUGAAAUUCCCGACCGGUGGUAUAGUCCAC
GAAAGUAUUUGCUUUGAUUUGGUGAAAUUCCAAAACCGACAGUAGAGUCUGGAUGAGAGAAGAUUC")

# Check number of high-scoring hits

length(searchHits)

# Extract the Rfam family accession and ID for the first hit

searchHits[[1]]$rfamAccession
searchHits[[1]]$rfamID
}
