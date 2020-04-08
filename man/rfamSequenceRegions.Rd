\name{rfamSequenceRegions}
\alias{rfamSequenceRegions}
\title{Retrieve all sequence regions belonging to an Rfam family
}
\description{
Retrieves all sequence regions encoding an RNA assigned to be a member of the
specified Rfam family. If a filename is provided, the list of sequence regions
will be saved to the path specified through filename as a tab-delimited file.
}
\usage{
rfamSequenceRegions(rfamFamily, filename=NULL)
}
\arguments{
  \item{rfamFamily}{string with an Rfam family accession or ID for which the
  member sequence regions should be retrieved.
  }
  \item{filename}{string indicating a path to which the retrieved sequence
  regions should be written as a tab-delimited file. In the default behaviour,
  filename=NULL, and no file will be written. Instead, the sequence regions will
  be returned as a data frame.
  }
}
\value{
A data frame with the sequence regions belonging to the specified Rfam family.
The data frame contains the following columns:
  \item{Sequence GenBank accession }{GenBank accession of the sequence with the
  sequence region belonging to the Rfam family}
  \item{Bit score }{Bit score for the alignment of the sequence region to the
  Rfam family}
  \item{Region start position }{Position of the sequence specified by the
  GenBank accession where the sequence region belonging to the Rfam family
  starts}
  \item{Region end position }{Position of the sequence specified by the
  GenBank accession where the sequence region belonging to the Rfam family
  ends}
  \item{Sequence description }{Description of the sequence where the sequence
  region belonging to the Rfam family is found}
  \item{Species }{Name of the species where the sequence region is present}
  \item{NCBI tax ID }{NCBI taxonomy ID of the species where the sequence region
  is present}
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
# Generate a data frame with the sequence regions belonging to the Rfam family
# with ID "FMN"

head(rfamSequenceRegions("FMN"))

# Generate a data frame with the sequence regions belonging to the Rfam family
# with accession "RF00360", and save the sequence regions as a tab-delimited
# text file.

head(rfamSequenceRegions("RF00360", filename="RF00360sequenceRegions.txt"))
\dontrun{
# Some Rfam families have too many associated sequence regions, in which case
# the list cannot be retrieved from the server.

rfamSequenceRegions("RF00174")
}
}
