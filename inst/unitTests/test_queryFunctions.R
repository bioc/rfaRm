library(RUnit)
library(rsvg)

## Test rfamFamilyAccessionToID

checkEquals(rfamFamilyAccessionToID("RF00050"), "FMN")

## Test rfamFamilyIDToAccession

checkEquals(rfamFamilyIDToAccession("FMN"), "RF00050")

## Test rfamFamilySummary

checkEquals(rfamFamilySummary("FMN")$rfamID, "FMN")

## Test rfamConsensusSecondaryStructure

checkTrue(grepl("^[]\\(\\[<\\{\\)>\\}\\.ABCDabcd]+$", rfamConsensusSecondaryStructure("FMN", format="DB")[2]))

## Test rfamSecondaryStructureXMLSVG

checkTrue(is.array(rsvg(charToRaw(rfamSecondaryStructureXMLSVG("FMN")))))

## Test rfamSecondaryStructurePlot

checkTrue(class(rfamSecondaryStructurePlot("FMN")) == "magick-image")

## Test rfamSeedAlignment

checkTrue(class(rfamSeedAlignment("FMN")) == "RNAMultipleAlignment")

## Test rfamSeedTree

checkTrue(grepl("\\[", rfamSeedTree("FMN", filename="testNHX.nhx")) & grepl("_", rfamSeedTree("FMN", filename="testNHX.nhx")) & is.character(rfamSeedTree("FMN", filename="testNHX.nhx")))

## Test rfamSeedTreeImage

checkTrue(class(rfamSeedTreeImage("FMN")) == "magick-image")

## Test rfamSequenceRegions

checkTrue(is.data.frame(rfamSequenceRegions("FMN")) & nrow(rfamSequenceRegions("FMN")) >= 1)

## Test rfamPDBMapping

checkTrue(is.data.frame(rfamPDBMapping("FMN")) & nrow(rfamPDBMapping("FMN")) >= 1)
