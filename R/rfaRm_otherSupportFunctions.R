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
    if(localOS == "Linux") {
        httrConfig <- config(ssl_cipher_list="DEFAULT@SECLEVEL=1")
        response <- with_config(config=httrConfig, POST(rfamApiSequenceSearchURL, 
                                                        accept_json(), body=queryBody,
                                                        encode="multipart"), user_agent("Safari/601.3.9"))
    }
    else {
        response <- POST(rfamApiSequenceSearchURL, accept_json(), body=queryBody,
                         encode="multipart", user_agent("Safari/601.3.9"))
    }
    return(content(response))
}

rfamSendSequenceSearchQuery2 <- function(sequence) {
    message("Running sequence search query. This might take a long time.")
    queryBody <- toJSON(list(
        query=sequence,
        databases=list("Rfam"),
        priority="high"
    ), auto_unbox = TRUE)
    if(localOS == "Linux") {
        httrConfig <- config(ssl_cipher_list="DEFAULT@SECLEVEL=1")
        response <- with_config(config=httrConfig, POST(rfamApiSequenceSearchURL2, 
                                                        accept_json(), body=queryBody,
                                                        encode="multipart"))
    }
    else {
        response <- POST(rfamApiSequenceSearchURL2, accept_json(), body=queryBody,
                         content_type_json())
    }
    return(content(response))
}

## Function to check the status of a previously sent query for a sequence-based
## search of the Rfam database.

rfamCheckSequenceSearchQuery <- function(responseURL) {
    if(localOS == "Linux") {
        httrConfig <- config(ssl_cipher_list="DEFAULT@SECLEVEL=1")
        response <- with_config(config=httrConfig, GET(responseURL, accept_json()))
    }
    else {
        response <- GET(responseURL, accept_json())
    }
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

rfamCheckSequenceSearchQuery2 <- function(responseURL) {
    if(localOS == "Linux") {
        httrConfig <- config(ssl_cipher_list="DEFAULT@SECLEVEL=1")
        response <- with_config(config=httrConfig, GET(responseURL, accept_json()))
    }
    else {
        response <- GET(responseURL, accept_json())
    }
    queryStatus <- content(response)$status
    return(queryStatus)
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
    if(localOS == "Linux") {
        httrConfig <- config(ssl_cipher_list="DEFAULT@SECLEVEL=1")
        result <- with_config(config=httrConfig, GET(responseURL, accept_json()))
    }
    else {
        result <- GET(responseURL, accept_json())
    }
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
                                                                                   "alignmentEndPositionHitSequence"=unlist(strsplit(hit$alignment$hit_seq, split=" +"))[length(unlist(strsplit(hit$alignment$hit_seq, split=" +")))],
                                                                                   "alignmentQuerySequence"=unlist(strsplit(hit$alignment$user_seq, split=" +"))[3],
                                                                                   "alignmentMatch"=trimws(substring(hit$alignment$match, 10), which="both"),
                                                                                   "alignmentHitSequence"=unlist(strsplit(hit$alignment$hit_seq, split=" +"))[3],
                                                                                   "alignmentSecondaryStructure"=unlist(strsplit(hit$alignment$ss, split=" +"))[2]))
    return(hitsList)
}

rfamRetrieveSequenceSearchResult2 <- function(resultURL) {
    if(localOS == "Linux") {
        httrConfig <- config(ssl_cipher_list="DEFAULT@SECLEVEL=1")
        result <- with_config(config=httrConfig, GET(resultURL, accept_json()))
    }
    else {
        result <- GET(resultURL, accept_json())
    }
    queryStatus <- status_code(result)
    if (queryStatus != 200) {
        stop("Sequence search has not completed")
    }
    autoparsedHits <- content(result)
    hitsList <- lapply(autoparsedHits, function(hit) {
        alignment <- hit$alignment
        splitAlignment <- strsplit(alignment, "\n")[[1]]
        secondaryStructure <- strsplit(splitAlignment[3], " +")[[1]][2]
        hitSequence <- strsplit(splitAlignment[4], " +")[[1]][4]
        querySequence <- strsplit(splitAlignment[6], " +")[[1]][4]
        alignmentMatch <- substring(splitAlignment[5], 15)
        list("rfamAccession"=hit$accession_rfam,
             "bitScore"=hit$score,
             "eValue"=hit$e_value,
             "strand"=hit$strand,
             "alignmentStartPositionQuerySequence"=hit$seq_from,
             "alignmentEndPositionQuerySequence"=hit$seq_to,
             "alignmentStartPositionHitSequence"=hit$mdl_from,
             "alignmentEndPositionHitSequence"=hit$mdl_to,
             "alignmentQuerySequence"=querySequence,
             "alignmentMatch"=alignmentMatch,
             "alignmentHitSequence"=hitSequence,
             "alignmentSecondaryStructure"=secondaryStructure)
        }
    )
    namesHits <- lapply(autoparsedHits, function(hit) {hit$target_name})
    names(hitsList) <- namesHits
    return(hitsList)
}

## Function to retrieve the complete set of clans defined in the Rfam database
## and the families belonging to each clan.

rfamGetClanDefinitions_old <- function() {
    clanHTMLTable <- html_table(xml_find_all(read_html(rfamClansListURL), "//table[@id]"))
    clanAccessions <- clanHTMLTable[[1]][,3][[1]]
    clanDefinitions <- vector(mode="list", length=length(clanAccessions))
    names(clanDefinitions) <- clanAccessions
    for (clan in clanAccessions) {
        clanFamiliesHTMLnodes <- xml_find_all(read_html(paste(rfamClanLookUpURL, clan, sep="")),
                                              "//span[@class='listItem']")
        clanFamilies <- regmatches(clanFamiliesHTMLnodes, regexpr("RF[0123456789]{5}", clanFamiliesHTMLnodes))
        clanDefinitions[[clan]] <- clanFamilies
    }
    return(clanDefinitions)
}

rfamGetClanDefinitions <- function() {
  clanMembershipCon <- gzcon(url("https://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/database_files/clan_membership.txt.gz"))
  clanMembershipTable <- read.delim(textConnection(readLines(clanMembershipCon)), header = FALSE)
  clanAccessions <- unique(clanMembershipTable[,1])
  clanDefinitions <- vector(mode="list", length=length(clanAccessions))
  names(clanDefinitions) <- clanAccessions
  for (clan in clanAccessions) {
    clanFamilies <- clanMembershipTable[clanMembershipTable[,1] == clan,2]
    clanDefinitions[[clan]] <- clanFamilies
  }
  return(clanDefinitions)
}


## Function to identify hits of a sequence search that overlap.

rfamFindOverlappingHits <- function(sequenceSearchHits, threshold) {
    querySequenceHitStarts <- as.integer(sapply(sequenceSearchHits, `[[`, "alignmentStartPositionQuerySequence"))
    querySequenceHitEnds <- as.integer(sapply(sequenceSearchHits, `[[`, "alignmentEndPositionQuerySequence"))
    for (startPoint in seq_len(length(querySequenceHitStarts))) {
        if (querySequenceHitStarts[startPoint] > querySequenceHitEnds[startPoint]) {
            oldEndPoint <- querySequenceHitEnds[startPoint]
            querySequenceHitEnds[startPoint] <- querySequenceHitStarts[startPoint]
            querySequenceHitStarts[startPoint] <- oldEndPoint
        }
    }
    querySequenceHitRanges <- IRanges(start=querySequenceHitStarts, end=querySequenceHitEnds)
    hitOverlaps <- findOverlaps(querySequenceHitRanges, drop.self=TRUE, drop.redundant=TRUE)
    overlapRegions <- pintersect(querySequenceHitRanges[queryHits(hitOverlaps)],
                                 querySequenceHitRanges[subjectHits(hitOverlaps)])
    overlapFractions <- width(overlapRegions)/pmin(width(querySequenceHitRanges[queryHits(hitOverlaps)]),
                                                   width(querySequenceHitRanges[subjectHits(hitOverlaps)]))
    overlappedHits <- as.matrix(hitOverlaps[overlapFractions >= threshold])
    colnames(overlappedHits) <- c("Hit1", "Hit2")
    return(overlappedHits)
}

## Function to remove overlapping hits of a sequence search that correspond
## to Rfam families belonging to the same clan.

rfamClanCompetitionFilter <- function(sequenceSearchHits, overlappedHits) {
    keepHit <- rep(TRUE, length(sequenceSearchHits))
    names(keepHit) <- names(sequenceSearchHits)
    for (overlap in seq_len(nrow(overlappedHits))) {
        hit1 <- sequenceSearchHits[[overlappedHits[overlap, 1]]]
        hit2 <- sequenceSearchHits[[overlappedHits[overlap, 2]]]
        clanHit1 <- names(rfamClanDefinitions[grep(hit1$rfamAccession, rfamClanDefinitions)])
        clanHit2 <- names(rfamClanDefinitions[grep(hit2$rfamAccession, rfamClanDefinitions)])
        if (isTRUE(hit1$rfamAccession == hit2$rfamAccession || clanHit1 == clanHit2)) {
            eValueHit1 <- hit1$eValue
            eValueHit2 <- hit2$eValue
            worseHit <- overlappedHits[overlap, which.max(c(eValueHit1, eValueHit2))]
            keepHit[worseHit] <- FALSE
        }
    }
    filteredHits <- sequenceSearchHits[keepHit]
    return(filteredHits)
}

## Function to check the query did not return an empty response.

checkEmptyResponse <- function(response) {
    if (length(content(response)) == 1) {
        if (is.na(content(response))) {
            stop("Requested data not currently available in the Rfam database")
        }
    }
}

