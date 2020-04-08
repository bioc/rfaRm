\name{rfamConsensusSecondaryStructure}
\alias{rfamConsensusSecondaryStructure}
\title{Retrieve the secondary structure of an Rfam family
}
\description{
Retrieves the consensus secondary structure and sequence of the specified Rfam
family. The notation of the secondary structure can be specified through the
"format" argument. If a filename is provided, the consensus secondary structure
and sequence will be saved in the specified format.
}
\usage{
rfamConsensusSecondaryStructure(rfamFamily, filename=NULL, format="DB")
}
\arguments{
  \item{rfamFamily}{string with an Rfam family accession or ID for which the
  consensus secondary structure and sequence should be retrieved.
  }
  \item{filename}{string indicating a path to which the retrieved consensus
  secondary structure and sequence should be written. The notation used for the
  RNA secondary structure is determined by the "format" argument. In the default
  behaviour, filename=NULL, and no file will be written. Instead, the consensus
  secondary structure and sequence will just be returned as a character vector.
  }
  \item{format}{string indicating the notation to be used for the RNA secondary
  structure. It can be either "DB" (extended Dot-Bracket notation; default) or
  "WUSS" (Washington University Secondary Structure notation).For a description
  of the different notations of RNA secondary structure, see \cr
  https://www.tbi.univie.ac.at/RNA/ViennaRNA/doc/html/rna_structure_notations.html.
  }
}
\value{
A character vector where the first element is the consensus sequence of the Rfam
family, and the second element is the consensus RNA secondary structure in the
specified format.
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
# Generate a character vector with the consensus secondary structure and
# sequence for the Rfam family with ID "FMN", with extended Dot-Bracket
# notation.

rfamConsensusSecondaryStructure("FMN", format="DB")

# Generate a character vector with the consensus secondary structure and
# sequence for the Rfam family with accession "RF00174", with extended
# Dot-Bracket notation, and save it to a file.

rfamConsensusSecondaryStructure("RF00174", filename="RF00174consensusSS.txt",
format="DB")

# Generate a character vector with the consensus secondary structure and
# sequence for the Rfam family with accession "RF00050", with WUSS notation, and
# save it to a file.

rfamConsensusSecondaryStructure("RF00050", filename="RF00050consensusSS.txt",
format="WUSS")
}
