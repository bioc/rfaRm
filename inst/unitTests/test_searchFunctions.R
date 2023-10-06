library(RUnit)

## Test rfamTextSearchFamilyAccession

checkTrue(is.element("RF00050", rfamTextSearchFamilyAccession("FMN")))

## Test rfamSequenceSearch

checkEquals(rfamSequenceSearch("GGAUCUUCGGGGCAGGGUGAAAUUCCCGACCGGUGGUAUAGUCCACGAAAGUAUUUGCUUUGAUUUGGUGAAAUUCCAAAACCGACAGUAGAGUCUGGAUGAGAGAAGAUUC")[[1]]$rfamAccession, "RF00050")

