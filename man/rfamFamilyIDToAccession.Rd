\name{rfamFamilyIDToAccession}
\alias{rfamFamilyIDToAccession}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Convert an Rfam family ID to an Rfam family accession}
\description{
Gets the Rfam family accession corresponding to the provided Rfam family ID.
}
\usage{
rfamFamilyIDToAccession(rfamFamilyID)
}
\arguments{
  \item{rfamFamilyID}{string with the Rfam family ID to be converted to an Rfam
  family accession.
}
}
\value{
A string with the Rfam family accession corresponding to the provided Rfam
family ID.
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
# Get the Rfam family accession corresponding to the Rfam family ID "FMN""

rfamFamilyIDToAccession("FMN")
}
