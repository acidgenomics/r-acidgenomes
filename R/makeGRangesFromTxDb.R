#' Make GRanges from TxDb object
#'
#' @note Updated 2021-01-22.
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
    keys <- AnnotationDbi::columns(object)
    switch(
        EXPR = level,
        "genes" = {
            what <- GenomicFeatures::genes
            columns <- "gene_id"
            args <- append(
                x = args,
                values = list(
                    "single.strand.genes.only" = TRUE
                )
            )
        },
        "transcripts" = {
            what <- GenomicFeatures::transcripts
            columns <- c("tx_id", "tx_name", "gene_id")
        }
    )
    args[["columns"]] <- columns
    suppressMessages({
        gr <- do.call(what = what, args = args)
    })
    assert(is(gr, "GRanges"))
    gr
}
