\name{rfamPDBMapping}
\alias{rfamPDBMapping}
\title{Retrieve PDB entries with the 3D structure of members of an Rfam family
}
\description{
Retrieves entries of the PDB database with the experimentally solved 3D
structure of members of the specified Rfam family, with correspondances between
residues of the PDB structure and positions in the covariance model of the Rfam
family. If a filename is provided, the list of PDB entries will be saved to the
path specified through filename as a tab-delimited file.
}
\usage{
rfamPDBMapping(rfamFamily, filename=NULL)
}
\arguments{
  \item{rfamFamily}{string with an Rfam family accession or ID for which the
  matching PDB entries should be retrieved.
  }
  \item{filename}{string indicating a path to which the retrieved PDB entries
  should be written as a tab-delimited file. In the default behaviour,
  filename=NULL, and no file will be written. Instead, the PDB entries will just
  be returned as a data frame.
  }
}
\value{
A data frame with PDB entries with experimentally solved 3D structures of
members of the specified Rfam family, as well as correspondances between
residues of the PDB structures and positions of the covariance model of the
Rfam family. The data frame contains the following columns:
  \item{Rfam accession }{Rfam accession of the Rfam family}
  \item{PDB ID }{ID of the PDB entry with a 3D structure of a member of the Rfam
  family}
  \item{Chain }{ID of the chain of the PDB entry that contains the coordinates
  for the 3D structure of a member of the Rfam family}
  \item{PDB start residue }{First residue of the 3D structure that can be
  matched to the sequence of a member of the Rfam family}
  \item{PDB end residue }{Last residue of the 3D structure that can be
  matched to the sequence of a member of the Rfam family}
  \item{CM start position }{First position of the covariance model of the Rfam
  family that can be matched to the 3D structure}
  \item{CM end position }{Last position of the covariance model of the Rfam
  family that can be matched to the 3D structure}
  \item{eValue }{Expectation value for the alignment of the sequence of the 3D
  structure to the Rfam family}
  \item{Bit score }{Bit score for the alignment of the sequence of the 3D
  structure to the Rfam family}
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
# Generate a data frame with the PDB structures of members of the Rfam family
# with ID "FMN"

rfamPDBMapping("FMN")

# Generate a data frame with the PDB structures of members of the Rfam family
# with accession "RF00174", and save the resulting PDB entries as a
# tab-delimited text file.

rfamPDBMapping("RF00174", filename="RF00174PDBentries.txt")
}
