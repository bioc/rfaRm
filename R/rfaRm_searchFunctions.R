## Function to search the Rfam database by keyword and retrieve the accessions
## of matching Rfam families

rfamTextSearchFamilyAccession <- function(query) {
    checkMultipleQuery(query)
    if (!is.character(query)) {
        stop("Input query is not a string")
    }
    result <- GET(rfamEbiApiURL,
                  query=list(query=query,
                             format="idlist"))
    matchList <- unlist(stri_split_lines(content(result, as="text")))
    return(matchList[grep("^RF", matchList)])
}

## Function to search the Rfam database by sequence and retrieve high-scoring
## hits with the consensus sequence of Rfam families. The sequence search
## is divided into 3 steps: sending the search query, checking for the search
## status and retrieving the results. Results are only retrieved after the
## search is confirmed to be finished. This is due to potentially long waiting
## times on the server side, which could lead to time-out errors.

rfamSequenceSearch <- function(sequence) {
    checkMultipleQuery(sequence)
    checkRNAString(sequence)
    sendQueryResult <- rfamSendSequenceSearchQuery(sequence)
    searchFinished <- FALSE
    while (!searchFinished) {
        Sys.sleep(10)
        checkQueryStatus <- rfamCheckSequenceSearchQuery(sendQueryResult$resultURL)
        if (checkQueryStatus == 200) {
            searchFinished <- TRUE
        }
    }
    searchResult <- rfamRetrieveSequenceSearchResult(sendQueryResult$resultURL)
    return(searchResult)
}
