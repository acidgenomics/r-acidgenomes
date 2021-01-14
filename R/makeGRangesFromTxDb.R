## FIXME NEED TO HANDLE THIS
##> genes <- genes(db)
##2714 genes were dropped because
##they have exons located on both
##strands of the same reference
##sequence or on more than one
##reference sequence, so cannot be
##represented by a single genomic
##range.
##Use
##'single.strand.genes.only=FALSE' to
##get all the genes in a GRangesList
##object, or use suppressMessages()
##to suppress this message.

## transcripts don't use an identifier, which is confusing.
## names should be defined as tx_name in this case instead.

## FIXME tx_name contains duplicates, so need to figure out how to handle.



## FIXME NEED TO RETHINK SUPPORT FOR REFSEQ HERE.

#' Make GRanges from TxDb object
#'
#' @note Updated 2021-01-14.
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
            args <- append(
                x = args,
                values = list(
                    "single.strand.genes.only" = TRUE
                )
            )
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
