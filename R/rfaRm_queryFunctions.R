## Function to convert an Rfam family accession to the corresponding Rfam
## family ID.

rfamFamilyAccessionToID <- function(rfamFamilyAccession) {
    checkMultipleQuery(rfamFamilyAccession)
    if (!grepl("^RF", rfamFamilyAccession)) {
        stop("Input is not a valid rfam accession")
    }
    if(localOS == "Linux") {
        httrConfig <- config(ssl_cipher_list="DEFAULT@SECLEVEL=1")
        result <- with_config(config=httrConfig, GET(paste(rfamApiBaseURL, rfamFamilyAccession, "/id", sep="" )))
    }
    else {
        result <- GET(paste(rfamApiBaseURL, rfamFamilyAccession, "/id", sep="" ))
    }
    checkEmptyResponse(result)
    return(content(result, as="text"))
}

## Function to convert an Rfam family ID to the corresponding Rfam family
## accession.

rfamFamilyIDToAccession <- function(rfamFamilyID) {
    checkMultipleQuery(rfamFamilyID)
    if(localOS == "Linux") {
        httrConfig <- config(ssl_cipher_list="DEFAULT@SECLEVEL=1")
        result <- with_config(config=httrConfig, GET(paste(rfamApiBaseURL, rfamFamilyID, "/acc", sep="" )))
    }
    else {
        result <- GET(paste(rfamApiBaseURL, rfamFamilyID, "/acc", sep="" ))
    }
    if(content(result, as="text") == "No such family") {
        stop("Input is not a valid rfam ID")
    }
    return(content(result, as="text"))
}

## Function to obtain a summary description of an Rfam family from its family
## accession or ID.

rfamFamilySummary <- function(rfamFamily) {
    checkMultipleQuery(rfamFamily)
    checkRfamEntry(rfamFamily)
    if(localOS == "Linux") {
        httrConfig <- config(ssl_cipher_list="DEFAULT@SECLEVEL=1")
        result <- with_config(config=httrConfig, GET(paste(rfamApiBaseURL, rfamFamily, sep="" ), accept_json()))
    }
    else {
        result <- GET(paste(rfamApiBaseURL, rfamFamily, sep="" ), accept_json())
    }
    checkEmptyResponse(result)
    autoparsedResult <- content(result)
    outputResult <- list("rfamReleaseNumber" = autoparsedResult$rfam$release$number,
                         "numberSequencesSeedAlignment" = autoparsedResult$rfam$curation$num_seed,
                         "sourceSeedAlignment" = autoparsedResult$rfam$curation$seed_source,
                         "numberSpecies" = autoparsedResult$rfam$curation$num_species,
                         "RNAType" = autoparsedResult$rfam$curation$type,
                         "numberSequencesAll" = autoparsedResult$rfam$curation$num_full,
                         "structureSource" = autoparsedResult$rfam$curation$structure_source,
                         "description" = autoparsedResult$rfam$description,
                         "rfamAccession" = autoparsedResult$rfam$acc,
                         "rfamID" = autoparsedResult$rfam$id,
                         "comment"= autoparsedResult$rfam$comment)
    return(outputResult)
}

## Function to obtain a representation of the secondary structure of an Rfam
## family in SVG format.

rfamSecondaryStructureXMLSVG <- function(rfamFamily, filename=NULL, plotType="norm") {
    checkMultipleQuery(rfamFamily)
    checkRfamEntry(rfamFamily)
    checkPlotType(plotType)
    if(localOS == "Linux") {
        httrConfig <- config(ssl_cipher_list="DEFAULT@SECLEVEL=1")
        result <- with_config(config=httrConfig, GET(paste(rfamApiBaseURL, rfamFamily, "/image/", plotType, sep=""), accept_xml()))
    }
    else {
        result <- GET(paste(rfamApiBaseURL, rfamFamily, "/image/", plotType, sep=""), accept_xml())
    }
    checkEmptyResponse(result)
    svg <- content(result, as="text")
    if (is.character(filename)){
        writeLines(svg, con=filename)
    }
    return(svg)
}

## Function to plot the secondary structure of an Rfam family in SVG format.
## Different plot types are available, listed in the documentation.

rfamSecondaryStructurePlot <- function(rfamFamily, filename=NULL, plotType="norm") {
    checkMultipleQuery(rfamFamily)
    checkRfamEntry(rfamFamily)
    checkPlotType(plotType)
    if(localOS == "Linux") {
        httrConfig <- config(ssl_cipher_list="DEFAULT@SECLEVEL=1")
        result <- with_config(config=httrConfig, GET(paste(rfamApiBaseURL, rfamFamily, "/image/", plotType, sep="")))
    }
    else {
        result <- GET(paste(rfamApiBaseURL, rfamFamily, "/image/", plotType, sep=""))
    }
    checkEmptyResponse(result)
    svg <- content(result)
    if (is.element(plotType, c("cons", "fcbp", "cov", "ent", "maxcm", "norm", "rchie"))) {
        image <- image_read(svg)
    }
    else if (is.element(plotType, c("rscape", "rscape-cyk"))) {
        image <- image_read(rsvg(svg))
    }
    if (is.null(filename)) {
        print(image)
    }
    else if (is.character(filename)){
        image_write(image, path=filename, format=unlist(strsplit(filename, "[.]"))[2])
    }
    return(image)
}

## Function to obtain the covariance model generated with Infernal of an Rfam
## family.

rfamCovarianceModel <- function(rfamFamily, filename=NULL) {
    checkMultipleQuery(rfamFamily)
    checkRfamEntry(rfamFamily)
    if(localOS == "Linux") {
        httrConfig <- config(ssl_cipher_list="DEFAULT@SECLEVEL=1")
        result <- with_config(config=httrConfig, GET(paste(rfamApiBaseURL, rfamFamily, "/cm", sep="")))
    }
    else {
        result <- GET(paste(rfamApiBaseURL, rfamFamily, "/cm", sep=""))
    }
    checkEmptyResponse(result)
    cmModel <- content(result, as="text")
    if (is.character(filename)) {
        writeLines(cmModel, con=filename)
    }
    return(cmModel)
}

## Function to obtain a list of sequence regions annotated as belonging to the
## specified Rfam family.

rfamSequenceRegions <- function(rfamFamily, filename=NULL) {
    checkMultipleQuery(rfamFamily)
    checkRfamEntry(rfamFamily)
    if(localOS == "Linux") {
        httrConfig <- config(ssl_cipher_list="DEFAULT@SECLEVEL=1")
        result <- with_config(config=httrConfig, GET(paste(rfamApiBaseURL, rfamFamily, 
                                                           "/regions", sep=""),
                                                     accept("text")))
    }
    else {
        result <- GET(paste(rfamApiBaseURL, rfamFamily, "/regions", sep=""),
                      accept("text"))
    }
    if (status_code(result) == 403) {
        stop("The family has too many regions to list.")
    }
    checkEmptyResponse(result)
    sequenceRegionsTable <- read.delim(text=content(result, as="text"), header=FALSE, skip=4)
    colnames(sequenceRegionsTable) <- c("Sequence GenBank accession",
                                        "Bit score",
                                        "Region start position",
                                        "Region end position",
                                        "Sequence description",
                                        "Species",
                                        "NCBI tax ID")
    if (is.character(filename)) {
        write.table(sequenceRegionsTable, file=file(filename), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
    }
    return(sequenceRegionsTable)
}

## Function to obtain the phylogenetic tree corresponding to the seed multiple
## sequence alignment used to define the specified Rfam family.

rfamSeedTree <- function(rfamFamily, filename) {
    checkMultipleQuery(rfamFamily)
    checkRfamEntry(rfamFamily)
    if(localOS == "Linux") {
        httrConfig <- config(ssl_cipher_list="DEFAULT@SECLEVEL=1")
        result <- with_config(config=httrConfig, GET(paste(rfamApiBaseURL, rfamFamily, 
                                                           "/tree", sep="")))
    }
    else {
        result <- GET(paste(rfamApiBaseURL, rfamFamily, "/tree", sep=""))
    }
    checkEmptyResponse(result)
    seedTree <- content(result, as="text")
    write(seedTree, file=filename)
    return(seedTree)
}

## Function to plot the phylogenetic tree corresponding to the seed multiple
## sequence alignment used to define the specified Rfam family.

rfamSeedTreeImage <- function(rfamFamily, filename=NULL, label="species") {
    checkMultipleQuery(rfamFamily)
    checkRfamEntry(rfamFamily)
    checkTreeLabel(label)
    if(localOS == "Linux") {
        httrConfig <- config(ssl_cipher_list="DEFAULT@SECLEVEL=1")
        result <- with_config(config=httrConfig, GET(paste(rfamApiBaseURL, rfamFamily, 
                                                           "/tree/label/", label, 
                                                           "/image", sep="")))
    }
    else {
        result <- GET(paste(rfamApiBaseURL, rfamFamily, "/tree/label/", label, "/image", sep=""))
    }
    checkEmptyResponse(result)
    gif <- content(result)
    image <- image_read(gif)
    if (is.null(filename)) {
        print(image)
    }
    else if (is.character(filename)){
        image_write(image, path=filename, format=unlist(strsplit(filename, "[.]"))[2])
    }
    return(image)
}

## Function to obtain a list of PDB entries with 3D structures corresponding to
## members of the specified Rfam family.

rfamPDBMapping <- function(rfamFamily, filename=NULL) {
    checkMultipleQuery(rfamFamily)
    checkRfamEntry(rfamFamily)
    if(localOS == "Linux") {
        httrConfig <- config(ssl_cipher_list="DEFAULT@SECLEVEL=1")
        result <- with_config(config=httrConfig, GET(paste(rfamApiBaseURL, rfamFamily, 
                                                           "/structures", sep=""),
                                                     accept_json()))
    }
    else {
        result <- GET(paste(rfamApiBaseURL, rfamFamily, "/structures", sep=""),
                      accept_json())
    }
    checkEmptyResponse(result)
    if (length(content(result)[[1]]) == 0){
        stop("No structures available in the PDB database for this family")
    }
    PDBTable <- rbindlist(content(result)[[1]])[, c(3, 6, 7, 9, 8, 5, 1, 4, 2)]
    colnames(PDBTable) <- c("Rfam accession",
                            "PDB ID",
                            "Chain",
                            "PDB start residue",
                            "PDB end residue",
                            "CM start position",
                            "CM end position",
                            "eValue",
                            "Bit score")
    if (is.character(filename)) {
        write.table(PDBTable, file=file(filename), sep="\t", quote=FALSE,
                    row.names=FALSE, col.names=TRUE)
    }
    return(PDBTable)
}

## Function to obtain the seed multiple sequence alignment used to defne the
## specified Rfam family.

rfamSeedAlignment <- function(rfamFamily, filename=NULL, format="stockholm") {
    checkMultipleQuery(rfamFamily)
    checkRfamEntry(rfamFamily)
    checkAlignmentFormat(format)
    if(localOS == "Linux") {
        httrConfig <- config(ssl_cipher_list="DEFAULT@SECLEVEL=1")
        result <- with_config(config=httrConfig, GET(paste(rfamApiBaseURL, rfamFamily, 
                                                           "/alignment/", format, sep="")))
    }
    else {
        result <- GET(paste(rfamApiBaseURL, rfamFamily, "/alignment/", format, sep=""))
    }
    checkEmptyResponse(result)
    alignment <- content(result)
    if (is.character(filename)) {
        writeLines(alignment, con=filename)
    }
    if (is.element(format, c("stockholm", "pfam"))) {
        return(readRNAMultipleAlignment(textConnection(alignment), "stockholm"))
    }
    else if (is.element(format, c("fasta", "fastau"))) {
        tmpFile <- tempfile()
        tmpFileConnection <- file(tmpFile)
        writeLines(alignment, con=tmpFileConnection)
        close(tmpFileConnection)
        return(readRNAMultipleAlignment(tmpFile, "fasta"))
    }
}

## Function to obtain the consensus secondary structure of the specified Rfam
## family.

rfamConsensusSecondaryStructure <- function(rfamFamily, filename=NULL, format="DB") {
    checkMultipleQuery(rfamFamily)
    checkRfamEntry(rfamFamily)
    if(localOS == "Linux") {
        httrConfig <- config(ssl_cipher_list="DEFAULT@SECLEVEL=1")
        result <- with_config(config=httrConfig, GET(paste(rfamApiBaseURL, rfamFamily, 
                                                           "/alignment", sep="")))
    }
    else {
        result <- GET(paste(rfamApiBaseURL, rfamFamily, "/alignment", sep=""))
    }
    checkEmptyResponse(result)
    alignment <- content(result)
    consensusSSLines <- grep("SS_cons", unlist(strsplit(alignment, "\n")), value=TRUE)
    consensusSSLines <- trimws(substring(consensusSSLines, 13))
    #consensusSS <- unlist(strsplit(consensusSSLine, " "))[length(unlist(strsplit(consensusSSLine, " ")))]
    #consensusSS <- unlist(strsplit(consensusSSLine, " "))[1]
    consensusSS <- paste(consensusSSLines, collapse="")
    if (format == "DB") {
        consensusSS <- WUSSToDB(consensusSS)
    }
    consensusSeqLines <- grep("GC RF", unlist(strsplit(alignment, "\n")), value=TRUE)
    consensusSeqLines <- trimws(substring(consensusSeqLines, 8))
    #consensusSeq <- unlist(strsplit(consensusSeqLine, " "))[length(unlist(strsplit(consensusSeqLine, " ")))]
    #consensusSeq <- unlist(strsplit(consensusSeqLine, " "))[1]
    consensusSeq <- paste(consensusSeqLines, collapse="")
    if (is.character(filename)) {
        writeLines(c(consensusSeq, consensusSS), con=filename)
    }
    return(c(consensusSeq, consensusSS))
}
