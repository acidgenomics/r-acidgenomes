## GENCODE
##
## Uses `gene_type` instead of `gene_biotype`.
## Note that `gene_id` and `gene_name` are nicely defined, so don't use `Name`.
## Consider removing gene and transcript versions automatically.
##
## GTF:
## >  [1] "source"                   "type"
## >  [3] "score"                    "phase"
## >  [5] "gene_id"                  "gene_type"
## >  [7] "gene_name"                "level"
## >  [9] "havana_gene"              "transcript_id"
## > [11] "transcript_type"          "transcript_name"
## > [13] "transcript_support_level" "tag"
## > [15] "havana_transcript"        "exon_number"
## > [17] "exon_id"                  "ont"
## > [19] "protein_id"               "ccdsid"
##
## GFF:
## >  [1] "source"                   "type"
## >  [3] "score"                    "phase"
## >  [5] "ID"                       "gene_id"
## >  [7] "gene_type"                "gene_name"
## >  [9] "level"                    "havana_gene"
## > [11] "Parent"                   "transcript_id"
## > [13] "transcript_type"          "transcript_name"
## > [15] "transcript_support_level" "tag"
## > [17] "havana_transcript"        "exon_number"
## > [19] "exon_id"                  "ont"
## > [21] "protein_id"               "ccdsid"



#' Detect PAR duplicates
#'
#' Match Ensembl spec by removing the duplicate PAR Y chromosome annotations.
#'
#' @note Updated 2020-01-20.
#' @noRd
.detectPARDupes <- function(object, idCol) {
    assert(is(object, "GRanges"))
    idCol <- match.arg(
        arg = idCol,
        choices = c("ID", "gene_id", "transcript_id")
    )
    dupes <- grep(pattern = "_PAR_Y$", x = mcols(object)[[idCol]], value = TRUE)
    if (hasLength(dupes)) {
        alertWarning(sprintf(
            "%d pseudoautosomal region (PAR) Y chromosome %s: {.var %s}.",
            length(dupes),
            ngettext(
                n = length(dupes),
                msg1 = "duplicate",
                msg2 = "duplicates"
            ),
            toString(dupes, width = 100L)
        ))
    }
    invisible(object)
}



.makeGenesFromGencodeGFF3 <- function(object) {
    assert(
        is(object, "GRanges"),
        isSubset(
            x = c("ID", "gene_biotype", "gene_id", "gene_name", "type"),
            y = colnames(mcols(object))
        )
    )
    object <- object[mcols(object)[["type"]] == "gene"]
    .detectPARDupes(object, idCol = "ID")
    object
}



.makeGenesFromGencodeGTF <- function(object) {
    object <- .makeGenesFromEnsemblGTF(object)
    .detectPARDupes(object, idCol = "gene_id")
    object
}



.makeTranscriptsFromGencodeGFF3 <- function(object) {
    assert(
        is(object, "GRanges"),
        isSubset(
            x = c(
                "gene_biotype", "gene_id",
                "transcript_biotype", "transcript_id"
            ),
            y = colnames(mcols(object))
        )
    )
    object <- object[mcols(object)[["type"]] == "transcript"]
    .detectPARDupes(object, idCol = "ID")
    object
}



.makeTranscriptsFromGencodeGTF <- function(object) {
    object <- .makeTranscriptsFromEnsemblGTF(object)
    .detectPARDupes(object, idCol = "transcript_id")
    object
}



.standardizeGencodeToEnsembl <- function(object) {
    assert(is(object, "GRanges"))
    mcolnames <- colnames(mcols(object))
    ## Match Ensembl spec, which uses `*_biotype` instead of `*_type`.
    mcolnames <- sub(
        pattern = "^gene_type$",
        replacement = "gene_biotype",
        x = mcolnames
    )
    mcolnames <- sub(
        pattern = "^transcript_type$",
        replacement = "transcript_biotype",
        x = mcolnames
    )
    colnames(mcols(object)) <- mcolnames
    object
}
