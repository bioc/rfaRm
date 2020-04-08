## Set of functions to check correctness of user input.

checkMultipleQuery <- function(inputString) {
    if (length(inputString) > 1) {
        stop("Multiple queries at a time not supported")
    }
}

checkRfamEntry <- function(inputString) {
    if (!grepl("^RF", inputString)) {
        checkID <- rfamFamilyIDToAccession(inputString)
        if(!grepl("^RF", checkID)) {
            stop("Please provide a valid Rfam accession or family ID")
        }
    }
    else {
        checkAccession <- rfamFamilyAccessionToID(inputString)
        if(checkAccession == "No such family") {
            stop("Please provide a valid Rfam accession or family ID")
        }
    }
}

checkPlotType <- function(plotType) {
    if(!is.element(plotType, c("cons", "fcbp", "cov", "ent", "maxcm", "norm", "rchie", "rscape", "rscape-cyk"))) {
        stop("Unrecognized plot type")
    }
}

checkTreeLabel <- function(treeLabelType) {
    if(!is.element(treeLabelType, c("species", "acc"))) {
        stop("Unrecognized tree labeling scheme. Please choose from species or accessions")
    }
}

checkAlignmentFormat <- function(alignmentFormat) {
    if(!is.element(alignmentFormat, c("stockholm", "pfam", "fasta", "fastau"))) {
        stop("Unrecognized or unsupported alignment format. Please choose from
             stockholm, pfam, fasta or fastau")
    }
}

checkSSformat <- function(SSFormat) {
    if(!is.element(SSFormat, c("DB", "WUSS"))) {
        stop("Unrecognized or unsupported secondary structure notation format")
    }
}

checkRNAString <- function(inputString) {
    if(!grepl('^[AUGCaugc]+$', gsub("[\r\n]", "", inputString)))
        stop("Please provide a string containing only standard RNA symbols")
}
