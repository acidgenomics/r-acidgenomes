## FIXME NEED TO RETHINK SUPPORT FOR REFSEQ HERE.

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
