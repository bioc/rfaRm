\name{rfamFamilySummary}
\alias{rfamFamilySummary}
\title{Get a summary describing an Rfam family
}
\description{
Gets a summary describing an Rfam family. The summary includes information
regarding the number of sequences and species contained in the family, a brief
description about the function of the family and the corresponding type of RNA.
}
\usage{
rfamFamilySummary(rfamFamily)
}
\arguments{
  \item{rfamFamily}{string with an Rfam family accession or ID for which a
  descriptive summary should be retrieved.
}
}
\value{
A list containing the following elements that describe the Rfam family:
  \item{rfamReleaseNumber }{Version of Rfam used in the query}
  \item{numberSequencesSeedAlignment }{Number of sequences used in the seed
  alignment of the Rfam family}
  \item{sourceSeedAlignment }{Published reference with the seed alignment of
  the Rfam family}
  \item{numberSpecies }{Number of species containing in their genomes sequences
  that belong to the Rfam family}
  \item{RNAType }{Keywords describing the type of RNA corresponding to the
  Rfam family}
  \item{numberSequencesAll }{Total number of sequences included in the Rfam
  family}
  \item{structureSource }{Description of the source for the secondary structure
  of the family (published, predicted) and PMID of the corresponding publication
  or prediction method if applicable}
  \item{description }{Descriptive name of the RNA included in the Rfam family }
  \item{rfamAccession }{Accession of the Rfam family }
  \item{rfamID }{ID of the Rfam family }
  \item{comment }{Short paragraph describing the characteristics and
  biological role of the Rfam family }
}
\references{
Ioanna Kalvari, Joanna Argasinska, Natalia Quinones-Olvera, Eric P Nawrocki,
Elena Rivas, Sean R Eddy, Alex Bateman, Robert D Finn, Anton I Petrov,
Rfam 13.0: shifting to a genome-centric resource for non-coding RNA families,
Nucleic Acids Research, Volume 46, Issue D1, 4 January 2018, Pages D335–D342,
https://doi.org/10.1093/nar/gkx1038

https://docs.rfam.org/en/latest/api.html
}
\examples{
# Get a summary for the Rfam family with ID "FMN""

rfamFamilySummary("FMN")

# Get a summary for the Rfam family with accession "RF00174""

rfamFamilySummary("RF00174")

# Get a brief paragraph describing the Rfam family with accession "RF00174""

rfamFamilySummary("RF00174")$comment
}

