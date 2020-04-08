\name{rfamSecondaryStructureXMLSVG}
\alias{rfamSecondaryStructureXMLSVG}
\title{Get an SVG file with a secondary structure representation
of an Rfam family
}
\description{
Gets a string corresponding to an SVG file (XML format) with a representation of
the secondary structure of the Rfam family, which can be read by rsvg or magick.
If a filename is provided, the SVG file will also be written to the specified
path.
}
\usage{
rfamSecondaryStructureXMLSVG(rfamFamily, filename=NULL, plotType="norm")
}
\arguments{
  \item{rfamFamily}{string with an Rfam family accession or ID for which a
  SVG file with a representation of the secondary structure should be retrieved.
  }
  \item{filename}{string indicating a path to which the SVG graphics should be
  written. In the default behaviour, filename=NULL, the SVG will not be written
  to any file, and instead will only be returned as a character string, which
  can be directly read by R packages such as rsvg and magick.
  }
  \item{plotType}{string indicating the desired type of secondary structure
  diagram. Can be one of "cons", "fcbp", "cov", "ent", "maxcm", "norm", "rchie",
  "rscape" or "rscape-cyk".  For a description of each type of diagram, see \cr
  https://docs.rfam.org/en/latest/api.html#secondary-structure-images.
  }
}
\value{
A string with the contents of a SVG file (XML format) with the representation
of a secondary structure diagram of the Rfam family. The string can be directly
read by rsvg or magick. If a filename is provided, the SVG file will also be
written to the specified path.
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
# Get a SVG representing a normal secondary structure diagram of the Rfam family
# with ID "FMN"

rfamSecondaryStructureXMLSVG("FMN")

# Get a SVG representing a basepair conservation secondary structure diagram
# for the Rfam family with accession "RF00174", and save it to a file

rfamSecondaryStructureXMLSVG("RF00174", filename="RF00174bpcons.svg", plotType="fcbp")
}
