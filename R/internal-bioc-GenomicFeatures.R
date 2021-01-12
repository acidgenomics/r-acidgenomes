#' Make TxDb from GRanges
#'
#' @details
#' This step can be noisy and generate expected warnings:
#'
#' ```
#' some exons are linked to transcripts not found in the file
#' The following orphan exon were dropped
#' The following orphan CDS were dropped
#' ```
#'
#' @seealso
#' - https://stackoverflow.com/questions/38603668
#' - https://stackoverflow.com/questions/16517795
#'
#' @noRd
.makeTxDbFromGFF <- function(file) {
    alert(sprintf(
        "Making {.var %s} from {.file %s} with {.pkg %s}::{.fun %s}.",
        "TxDb", file,
        "GenomicFeatures", "makeTxDbFromGFF"
    ))
    requireNamespaces("GenomicFeatures")
    assert(isAFile(file))
    txdb <- withCallingHandlers(expr = {
        tryCatch(
            expr = {
                GenomicFeatures::makeTxDbFromGFF(file)
            },
            error = function(e) {
                stop(paste(
                    "Failed to make TxDb using ",
                    "'GenomicFeatures::makeTxDbFromGFF()'",
                    sep = "\n"
                ))
            }
        )
    }, warning = function(w) {
        if (grepl(pattern = "stop_codon", x = conditionMessage(w))) {
            invokeRestart("muffleWarning")
        }
    })
    assert(is(txdb, "TxDb"))
    validObject(txdb)
    txdb
}



#' Make GRanges from TxDb object
#'
#' @note Updated 2021-01-12.
#' @noRd
#'
#' @return `GRanges`.
.makeGRangesFromTxDb <- function(
    object,
    level = c("genes", "transcripts")
) {
    assert(is(object, "TxDb"))
    level <- match.arg(level)
    args <- list("x" = object)
    ## FIXME THIS RETURNS INCORRECTLY FOR DROSOPHILA.
    ## FIXME NEED TO PARSE THE GTF DIRECTLY TO GET METADATA...
    keys <- AnnotationDbi::columns(object)
    switch(
        EXPR = level,
        "genes" = {
            fun <- GenomicFeatures::genes
            columns <- "gene_id"
        },
        "transcripts" = {
            fun <- GenomicFeatures::transcripts
            columns <- c("tx_id", "tx_name", "gene_id")
        }
    )
    args[["columns"]] <- columns
    gr <- do.call(what = fun, args = args)
    assert(is(gr, "GRanges"))
    gr
}
