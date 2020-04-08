\name{rfamCovarianceModel}
\alias{rfamCovarianceModel}
\title{Get the covariance model associated to an Rfam family
}
\description{
Gets the covariance model generated with the Infernal software of the Rfam
family. Covariance models can be used in Infernal to search through a list of
sequences for high-scoring hits to the models. If a filename is provided, the
covariance model will be saved to the specified path.
}
\usage{
rfamCovarianceModel(rfamFamily, filename=NULL)
}
\arguments{
  \item{rfamFamily}{string with an Rfam family accession or ID for which the
  covariance model should be retrieved.
  }
  \item{filename}{string indicating a path to which the covariance model should
  be written. In the default behaviour, filename=NULL, the covariance model will
  not be saved into any file, and instead will only be returned as a string.
  }
}
\value{
A string with the covariance model of the Rfam family in the Infernal format.
If a filename is provided, the covariance model file will also be written to the
specified path.
}
\references{
Ioanna Kalvari, Joanna Argasinska, Natalia Quinones-Olvera, Eric P Nawrocki,
Elena Rivas, Sean R Eddy, Alex Bateman, Robert D Finn, Anton I Petrov,
Rfam 13.0: shifting to a genome-centric resource for non-coding RNA families,
Nucleic Acids Research, Volume 46, Issue D1, 4 January 2018, Pages D335–D342,
https://doi.org/10.1093/nar/gkx1038


Eric P Nawrocki, Sean R Eddy, Infernal 1.1: 100-fold faster RNA homology
searches, Bioinformatics, Volume 29, Issue 22, 15 November 2013,
Pages 2933–2935, https://doi.org/10.1093/bioinformatics/btt509

https://docs.rfam.org/en/latest/api.html

http://eddylab.org/infernal/Userguide.pdf
}
\examples{
# Retrieve the covariance model for the Rfam family with ID "FMN" and save it
# to a file

rfamCovarianceModel("FMN", filename="FMNcovarianceModel.cm")

# Retrieve the covariance model for the Rfam family with accession "RF00174" and
# save it to a file

rfamCovarianceModel("RF00174", filename="RF00174covarianceModel.cm")
}

