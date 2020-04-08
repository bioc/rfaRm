\name{rfamTextSearchFamilyAccession}
\alias{rfamTextSearchFamilyAccession}
\title{Search the Rfam database by a given keyword}
\description{
Searches the Rfam database for entries containing a specified keyword, such as
"ribozyme" or "FMN", and returns the accession numbers of the Rfam families
containing the keyword in its family ID, summary or description.
}
\usage{
rfamTextSearchFamilyAccession(query)
}
\arguments{
  \item{query}{string with the keyword to be searched in the Rfam database.
  }
}
\value{
A character vector with the accessions of the Rfam families matching the search.
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
# Search Rfam families associated to the keyword "FMN"

rfamTextSearchFamilyAccession("FMN")
}
