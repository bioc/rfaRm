\name{rfamFamilyAccessionToID}
\alias{rfamFamilyAccessionToID}
\title{Convert an Rfam family accession to an Rfam family ID}
\description{
Gets the Rfam family ID corresponding to the provided Rfam family accession.
}
\usage{
rfamFamilyAccessionToID(rfamFamilyAccession)
}
\arguments{
  \item{rfamFamilyAccession}{string with the Rfam family accession to be converted to an Rfam
  family ID.
}
}
\value{
A string with the Rfam family ID corresponding to the provided Rfam family
accesion. The Rfam family ID usually describes the type of RNA associated to
the Rfam family.
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
# Get the Rfam family ID corresponding to the Rfam accession RF00050

rfamFamilyAccessionToID("RF00050")
}
