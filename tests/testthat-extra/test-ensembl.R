context("Ensembl annotations")

## FIXME Need to ensure these run.
## FIXME Need to test GRCm39 support.

map <- import(
    file = system.file(
        "extdata", "map-genome-build.rds",
        package = packageName()
    ),
    quiet = TRUE
)
keep <- as.logical(map[["ensembldb"]])
map <- map[keep, , drop = FALSE]

## Note that we're running at transcript level here to check the gene merge.
##
## Potential issues:
## - *Caenorhabditis elegans*
##     - Invalid transcript identifiers.
## - *Canis familiaris*
##     - Using the full *Canis lupus familiaris* won't match.
## - *Saccharomyces cerevisiae*
##     - Invalid gene identifiers.
##     - Iinvalid transcript identifiers.
##     - No gene names.

test_that("UCSC genome build remaps", {
    mapply(
        organism = map[["organism"]],
        genomeBuild = map[["ucsc"]],
        FUN = function(organism, genomeBuild) {
            object <- makeGRangesFromEnsembl(
                organism = organism,
                genomeBuild = genomeBuild,
                level = "transcripts"
            )
            expect_s4_class(object, "GRanges")
        },
        SIMPLIFY = FALSE
    )
})
