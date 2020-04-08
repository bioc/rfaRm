\name{rfamSecondaryStructurePlot}
\alias{rfamSecondaryStructurePlot}
\title{Generate a diagram of the secondary structure of an Rfam family
}
\description{
Plots a diagram of the specified type of the secondary structure of an Rfam
family. If a filename is provided, the diagram will be saved to the specified
path.
}
\usage{
rfamSecondaryStructurePlot(rfamFamily, filename=NULL, plotType="norm")
}
\arguments{
  \item{rfamFamily}{string with an Rfam family accession or ID for which the
  secondary structure should be plotted.
  }
  \item{filename}{string indicating a path to which the plotted diagram should
  be written. In the default behaviour, filename=NULL, the diagram will not be
  saved into any file, and instead will only be plotted on R's graphics display.
  Must be a file type supported by magick, such as PNG, JPEG, GIF, RGB or RGBA.
  }
  \item{plotType}{string indicating the desired type of secondary structure
  diagram. Can be one of "cons", "fcbp", "cov", "ent", "maxcm", "norm"
  (default), "rchie", "rscape" or "rscape-cyk".  For a description of each type
  of diagram, see \cr
  https://docs.rfam.org/en/latest/api.html#secondary-structure-images.
  }
}
\value{
A magick image object with the requested diagram. If a filename was specified,
the diagram plot will also be saved to the specified path.
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
# Plot a normal secondary structure diagram of the Rfam family with ID "FMN"

rfamSecondaryStructurePlot("FMN")

# Save a basepair conservation secondary structure diagram of the Rfam family
# with accession "RF00174" as a PNG file.

rfamSecondaryStructurePlot("RF00174", filename="RF00174bpcons.png", plotType="fcbp")
}
