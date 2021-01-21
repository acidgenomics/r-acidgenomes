## FlyBase
##
## GTF:
## > [1] "source"            "type"              "score"
## > [4] "phase"             "gene_id"           "gene_symbol"
## > [7] "transcript_id"     "transcript_symbol" "#"



## Compatible with Ensembl importer after we run `.standardizeFlyBaseGFF()`,
## which is called in `.makeGenesFromGFF()`.
.makeGenesFromFlyBaseGTF <- function(object) {
    assert(is(object, "GRanges"))
    object <- .makeGenesFromEnsemblGTF(object)
    object
}



.makeTranscriptsFromFlyBaseGTF <- function(object) {
    assert(is(object, "GRanges"))
    ## Note that FlyBase uses non-standard transcript types.
    keep <- grepl(
        pattern = paste(c("^pseudogene$", "RNA$"), collapse = "|"),
        x = mcols(object)[["type"]],
        ignore.case = TRUE
    )
    object <- object[keep]
    object
}



.standardizeFlyBaseToEnsembl <- function(object) {
    assert(is(object, "GRanges"))
    mcolnames <- colnames(mcols(object))
    ## Match Ensembl spec by renaming `*_symbol` to `*_name`.
    mcolnames <- sub(
        pattern = "^gene_symbol$",
        replacement = "gene_name",
        x = mcolnames
    )
    mcolnames <- sub(
        pattern = "^transcript_symbol$",
        replacement = "transcript_name",
        x = mcolnames
    )
    colnames(mcols(object)) <- mcolnames
    object
}
