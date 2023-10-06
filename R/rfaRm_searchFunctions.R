## Function to search the Rfam database by keyword and retrieve the accessions
## of matching Rfam families

rfamTextSearchFamilyAccession <- function(query) {
    checkMultipleQuery(query)
    if (!is.character(query)) {
        stop("Input query is not a string")
    }
    if(localOS == "Linux") {
        httrConfig <- config(ssl_cipher_list="DEFAULT@SECLEVEL=1")
        result <- with_config(config=httrConfig, GET(rfamEbiApiURL,
                                                     query=list(query=query,
                                                                format="idlist")))
    }
    else {
        result <- GET(rfamEbiApiURL,
                      query=list(query=query,
                                 format="idlist"))
    }
    matchList <- unlist(stri_split_lines(content(result, as="text")))
    return(matchList[grep("^RF", matchList)])
}

## Function to search the Rfam database by sequence and retrieve high-scoring
## hits with the consensus sequence of Rfam families. The sequence search
## is divided into 3 steps: sending the search query, checking for the search
## status and retrieving the results. Results are only retrieved after the
## search is confirmed to be finished. This is due to potentially long waiting
## times on the server side, which could lead to time-out errors.

rfamSequenceSearch_oldApi <- function(sequence, fragmentsOverlap=1000, clanCompetitionFilter=TRUE, clanOverlapThreshold=0.5) {
    checkMultipleQuery(sequence)
    checkRNAString(sequence)
    if (!(isTRUE(clanCompetitionFilter) | isFALSE(clanCompetitionFilter))) {
        stop("clanCompetitionFilter value must be either TRUE or FALSE.")
    }
    if (clanCompetitionFilter & !(clanOverlapThreshold >= 0 & clanOverlapThreshold <= 1)) {
        stop("clanOverlapThreshold should be a number between 0 and 1.")
    }
    if (nchar(sequence) > 10000) {
        if (!(fragmentsOverlap%%1 == 0) | fragmentsOverlap < 1) {
            stop("fragmentsOverlap should be a positive integer.")
        }
        fragmentEndPoints <- c(seq(from=10000, to=nchar(sequence), by=10000-fragmentsOverlap), nchar(sequence))
        fragmentStartPoints <- seq(from=1, by=10000-fragmentsOverlap, length.out=length(fragmentEndPoints))
        splitSequence <- character(length(fragmentEndPoints))
        splitSequence <- vapply(seq_len(length(fragmentEndPoints)),
                                FUN=function(x) substring(sequence, fragmentStartPoints[x], fragmentEndPoints[x]),
                                FUN.VALUE=character(1))
    }
    else {
        fragmentStartPoints <- 1
        fragmentEndPoints <- 10000
        splitSequence <- sequence
    }
    fullResults <- vector(mode="list", length=length(splitSequence))
    for (fragment in seq_len(length(splitSequence))) {
        if (length(splitSequence) > 1) {
            message(paste("Processing fragment", fragment, sep=" "))
        }
        sendQueryResult <- rfamSendSequenceSearchQuery(splitSequence[fragment])
        searchFinished <- FALSE
        while (!searchFinished) {
            Sys.sleep(10)
            checkQueryStatus <- rfamCheckSequenceSearchQuery(sendQueryResult$resultURL)
            if (checkQueryStatus == 200) {
                searchFinished <- TRUE
            }
        }
        searchResult <- rfamRetrieveSequenceSearchResult(sendQueryResult$resultURL)
        for (hit in seq_len(length(searchResult))) {
            searchResult[[hit]]$alignmentStartPositionQuerySequence <- as.character(as.integer(searchResult[[hit]]$alignmentStartPositionQuerySequence)+fragmentStartPoints[fragment])
            searchResult[[hit]]$alignmentEndPositionQuerySequence <- as.character(as.integer(searchResult[[hit]]$alignmentEndPositionQuerySequence)+fragmentStartPoints[fragment])
        }
        fullResults[[fragment]] <- searchResult
    }
    fullResults <- unlist(fullResults, recursive=FALSE)
    if (clanCompetitionFilter & (length(fullResults) >1)) {
        overlappingHits <- rfamFindOverlappingHits(fullResults, clanOverlapThreshold)
        fullResults <- rfamClanCompetitionFilter(fullResults, overlappingHits)
    }
    return(fullResults)
}


rfamSequenceSearch <- function(sequence, fragmentsOverlap=1000, clanCompetitionFilter=TRUE, clanOverlapThreshold=0.5) {
    checkMultipleQuery(sequence)
    checkRNAString(sequence)
    if (!(isTRUE(clanCompetitionFilter) | isFALSE(clanCompetitionFilter))) {
        stop("clanCompetitionFilter value must be either TRUE or FALSE.")
    }
    if (clanCompetitionFilter & !(clanOverlapThreshold >= 0 & clanOverlapThreshold <= 1)) {
        stop("clanOverlapThreshold should be a number between 0 and 1.")
    }
    if (nchar(sequence) > 10000) {
        if (!(fragmentsOverlap%%1 == 0) | fragmentsOverlap < 1) {
            stop("fragmentsOverlap should be a positive integer.")
        }
        fragmentEndPoints <- c(seq(from=10000, to=nchar(sequence), by=10000-fragmentsOverlap), nchar(sequence))
        fragmentStartPoints <- seq(from=1, by=10000-fragmentsOverlap, length.out=length(fragmentEndPoints))
        splitSequence <- character(length(fragmentEndPoints))
        splitSequence <- vapply(seq_len(length(fragmentEndPoints)),
                                FUN=function(x) substring(sequence, fragmentStartPoints[x], fragmentEndPoints[x]),
                                FUN.VALUE=character(1))
    }
    else {
        fragmentStartPoints <- 1
        fragmentEndPoints <- 10000
        splitSequence <- sequence
    }
    fullResults <- vector(mode="list", length=length(splitSequence))
    for (fragment in seq_len(length(splitSequence))) {
        if (length(splitSequence) > 1) {
            message(paste("Processing fragment", fragment, sep=" "))
        }
        sendQueryResult <- rfamSendSequenceSearchQuery2(splitSequence[fragment])
        searchFinished <- FALSE
        while (!searchFinished) {
            Sys.sleep(10)
            checkQueryStatus <- rfamCheckSequenceSearchQuery2(paste(rfamApiSequenceSearchURL2Check,
                                                                    sendQueryResult$job_id,
                                                                    sep=""))
            if (checkQueryStatus == "success") {
                searchFinished <- TRUE
            }
        }
        searchResult <- rfamRetrieveSequenceSearchResult2(paste(rfamApiSequenceSearchURL2Retrieve,
                                                                sendQueryResult$job_id,
                                                                sep=""))
        for (hit in seq_len(length(searchResult))) {
            searchResult[[hit]]$alignmentStartPositionQuerySequence <- as.character(as.integer(searchResult[[hit]]$alignmentStartPositionQuerySequence)+fragmentStartPoints[fragment])
            searchResult[[hit]]$alignmentEndPositionQuerySequence <- as.character(as.integer(searchResult[[hit]]$alignmentEndPositionQuerySequence)+fragmentStartPoints[fragment])
        }
        fullResults[[fragment]] <- searchResult
    }
    fullResults <- unlist(fullResults, recursive=FALSE)
    if (clanCompetitionFilter & (length(fullResults) >1)) {
        overlappingHits <- rfamFindOverlappingHits(fullResults, clanOverlapThreshold)
        fullResults <- rfamClanCompetitionFilter(fullResults, overlappingHits)
    }
    return(fullResults)
}

