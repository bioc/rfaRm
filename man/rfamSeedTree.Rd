\name{rfamSeedTree}
\alias{rfamSeedTree}
\title{Get the phylogenetic tree of the seed multiple alignment of an
Rfam family
}
\description{
Gets the phylogenetic tree of the seed multiple alignment associated to the
Rfam family. The tree is retrieved in the NHX format (New Hampshire extended
format. The tree will be saved to the path specified through filename. The
generated file can be read directly with the treeio R package.
}
\usage{
rfamSeedTree(rfamFamily, filename)
}
\arguments{
  \item{rfamFamily}{string with an Rfam family accession or ID for which the
  phylogenetic tree of the seed alignment should be retrieved.
  }
  \item{filename}{string indicating a path to which the phylogenetic tree in the
  NHX format should be written.
  }
}
\value{
A string with the phylogenetic tree of the seed alignment of the Rfam family in
the NHX format. The tree will also be written to the specified path, producing a
file that can be read by the treeio package.
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
# Retrieve the phylogenetic tree of the seed alignment for the Rfam family with
# ID "FMN" and save it to a file

rfamSeedTree("FMN", filename="FMNseedTree.nhx")

# Retrieve the phylogenetic tree of the seed alignment for the Rfam family with
# accession "RF00174" and save it to a file

rfamSeedTree("RF00174", filename="RF00174seedTree.nhx")
}
