## Function to convert RNA secondary structure string from the WUSS notation  to
## the extended Dot-Bracket notation.

WUSSToDB <- function(WUSSString) {
    WUSSVector <- unlist(strsplit(WUSSString, ""))
    DBVector <- character(length=length(WUSSVector))
    for (i in seq_len(length(WUSSVector))) {
        if (is.element(WUSSVector[i], c("(", "[", "<", "A", "B", "C", "D", ")",
                                        "]", ">", "a", "b", "c", "d"))) {
            DBVector[i] <- WUSSVector[i]
        }
        else if (is.element(WUSSVector[i], c("_", ":", ",", ".", "-"))) {
            DBVector[i] <- "."
        }
        else {
            stop("Invalid characters present in WUSS string")
        }
    }
    return(paste(DBVector, collapse=""))
}

## Function to send the query for a sequence-based search to the Rfam database.

rfamSendSequenceSearchQuery <- function(sequence) {
    message("Running sequence search query. This might take a long time.")
    queryBody <- list(sequence, "Submit")
    names(queryBody) <- c("seq", "submit")
    response <- POST(rfamApiSequenceSearchURL, accept_json(), body=queryBody,
                     encode="multipart")
    return(content(response))
}

## Function to check the status of a previously sent query for a sequence-based
## search of the Rfam database.

rfamCheckSequenceSearchQuery <- function(responseURL) {
    response <- GET(responseURL, accept_json())
    queryStatus <- status_code(response)
    if (queryStatus == 202) {
        message("Sequence search is running, please wait.")
        return(queryStatus)
    }
    else if (queryStatus == 200) {
        message("Sequence search completed successfully.")
        return(queryStatus)
    }
    else {
        message("Malformed query or server unavailable. Please try again.")
        return(queryStatus)
    }
}

## Function to retrieve the results of a completed query for a sequence-based
## search of the Rfam database.

rfamRetrieveSequenceSearchResult <- function(responseURL) {
    result <- GET(responseURL, accept_json())
    queryStatus <- status_code(result)
    if (queryStatus != 200) {
        stop("Sequence search has not completed")
    }
    autoparsedHits <- content(result)$hits
    hitsList <- lapply(unlist(autoparsedHits, recursive=FALSE), function(hit) list("rfamAccession"=hit$acc,
                                                                                   "rfamID"=hit$id,
                                                                                   "bitScore"=hit$score,
                                                                                   "eValue"=hit$E,
                                                                                   "strand"=hit$strand,
                                                                                   "alignmentStartPositionQuerySequence"=hit$start,
                                                                                   "alignmentEndPositionQuerySequence"=hit$end,
                                                                                   "alignmentStartPositionHitSequence"=unlist(strsplit(hit$alignment$hit_seq, split=" +"))[2],
                                                                                   "alignmentEndPositionHitSequence"=unlist(strsplit(hit$alignment$hit_seq, split=" +"))[4],
                                                                                   "alignmentQuerySequence"=unlist(strsplit(hit$alignment$user_seq, split=" +"))[3],
                                                                                   "alignmentMatch"=trimws(substring(hit$alignment$match, 10), which="both"),
                                                                                   "alignmentHitSequence"=unlist(strsplit(hit$alignment$hit_seq, split=" +"))[3],
                                                                                   "alignmentSecondaryStructure"=unlist(strsplit(hit$alignment$ss, split=" +"))[2]))
    return(hitsList)
}

## Function to check the query did not return an empty response.

checkEmptyResponse <- function(response) {
    if (length(content(response)) == 1) {
        if (is.na(content(response))) {
            stop("Requested data not currently available in the Rfam database")
        }
    }
}
