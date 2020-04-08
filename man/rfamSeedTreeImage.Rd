\name{rfamSeedTreeImage}
\alias{rfamSeedTreeImage}
\title{Get an image of the phylogenetic tree of the seed multiple alignment of
an Rfam family
}
\description{
Plots an image of the phylogenetic tree of the seed multiple alignment
associated to the Rfam family. If a filename is provided, the image will be
saved to the path specified through filename. The tree can be labeled with
species names or sequence accessions.
}
\usage{
rfamSeedTreeImage(rfamFamily, filename=NULL, label="species")
}
\arguments{
  \item{rfamFamily}{string with an Rfam family accession or ID for which the
  phylogenetic tree of the seed alignment should be plotted.
  }
  \item{filename}{string indicating a path to which the plot of the phylogenetic
  tree should be written.
  }
  \item{label}{string indicating the type of labels that should be added to
  the plot of the phylogenetic tree. Should be either "species" for labeling
  with species names (default) or "acc" for labeling with sequence accessions.
  }
}
\value{
A magick image object with the requested phylogenetic tree plot in GIF format.
If a filename was specified, the phylogenetic tree plot will also be saved to
the specified path.
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
# Display a plot of the phylogenetic tree of the seed alignment for the Rfam
# family with ID "FMN" labeled with species names

rfamSeedTreeImage("FMN", label="species")

# Display a plot of the phylogenetic tree of the seed alignment for the Rfam
# family with accession "RF00147" labeled with sequence accessions and save it
# to a PNG file.

rfamSeedTreeImage("RF00147", filename="RF00147seedTreeImage.png", label="acc")
}
