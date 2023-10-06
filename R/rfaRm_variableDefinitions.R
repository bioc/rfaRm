## URL for keyword-based searches of the Rfam database
rfamEbiApiURL <- 'https://www.ebi.ac.uk/ebisearch/ws/rest/rfam'

## Base URL used for most queries to the Rfam database
rfamApiBaseURL <- 'https://rfam.org/family/'

## Base URL used for sequence-based searches of the Rfam database
rfamApiSequenceSearchURL <- 'https://rfam.org/search/sequence'
rfamApiSequenceSearchURL2 <- 'https://search.rnacentral.org/api/submit-job'
rfamApiSequenceSearchURL2Check <- "https://search.rnacentral.org/api/infernal-status/"
rfamApiSequenceSearchURL2Retrieve <- "https://search.rnacentral.org/api/infernal-result/"

## URL used to retrieve the list of clans defined in the Rfam database
rfamClansListURL <- 'https://rfam.org/clans'

## Base URL used to retrieve the Rfam families belonging to a specific clan
rfamClanLookUpURL <- 'https://rfam.org/clan/'

## OS
localOS <- Sys.info()["sysname"]

## Named list where the name of each element is an Rfam clan, and each element
## is a character vector with the Rfam families of the corresponding clan
rfamClanDefinitions <- NULL
