\name{rfamSeedAlignment}
\alias{rfamSeedAlignment}
\title{Retrieve the seed alignment of an Rfam family
}
\description{
Retrieves the seed multiple alignment of the specified Rfam family. The seed
alignment is used to determine the covariance model defining each Rfam family,
and comprises only a subset of all the members of each family. If a filename is
provided, the alignment will be saved in the specified format.
}
\usage{
rfamSeedAlignment(rfamFamily, filename=NULL, format="stockholm")
}
\arguments{
  \item{rfamFamily}{string with an Rfam family accession or ID for which the
  seed alignment should be retrieved.
  }
  \item{filename}{string indicating a path to which the retrieved seed alignment
  should be written. The format of the written file is determined by the
  "format" argument. In the default behaviour, filename=NULL, and no file will
  be written. Instead, the seed alignment will just be returned as a Biostrings
  MultipleAlignment object.
  }
  \item{format}{string indicating the format into which the seed alignment
  should be retrieved and written to filename Can be one of "stockholm"
  (default), "pfam", "fasta" or "fastau". For a description of each type of
  format, see \cr
  https://docs.rfam.org/en/latest/api.html#stockholm-format-alignment.
  }
}
\value{
A Biostrings MultipleAlignment object with the aligned RNA sequences included in
the seed alignment. The names of the sequences contain their GenBank accession
numbers and the start and end positions of the region corresponding to the
sequence region member of the Rfam family.
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
# Generate a Biostrings MultipleAlignment object with the aligned RNA sequences
# used to generate the seed alignment for the Rfam family with ID "FMN"

rfamSeedAlignment("FMN")

# Generate a Biostrings MultipleAlignment object with the aligned RNA sequences
# used to generate the seed alignment for the Rfam family with accession
# "RF00174", and save it to a file in Stockholm format

rfamSeedAlignment("RF00174", filename="RF00174seedAlignment.stk", format="stockholm")

# Generate a Biostrings MultipleAlignment object with the aligned RNA sequences
# used to generate the seed alignment for the Rfam family with accession
# "RF00174", and save it to a file in FASTA format

rfamSeedAlignment("RF00174", filename="RF00174seedAlignment.fasta", format="fasta")
}
