#' Make GRanges from TxDb object
#'
#' @export
#' @note Updated 2021-01-25.
#'
#' @return `GRanges`.
#'
#' @examples
#' file <- pasteURL(
#'     "ftp.ensembl.org",
#'     "pub",
#'     "release-102",
#'     "gtf",
#'     "homo_sapiens",
#'     "Homo_sapiens.GRCh38.102.gtf.gz",
#'     protocol = "ftp"
#' )
#' txdb <- makeTxDbFromGFF(file)
#' gr <- makeGRangesFromTxDb(txdb)
#' print(gr)
makeGRangesFromTxDb <- function(
    object,
    level = c("transcripts", "genes", "exons", "cds"),
    ignoreVersion = TRUE
) {
    assert(is(object, "TxDb"))
    level <- match.arg(level)
    keys <- columns(object)
    colsList <- list(
        "cds" = grep(pattern = "^CDS", x = keys, value = TRUE),
        "exon" = grep(pattern = "^EXON", x = keys, value = TRUE),
        "gene" = grep(pattern = "^GENE", x = keys, value = TRUE),
        "tx" = grep(pattern = "^TX", x = keys, value = TRUE)
    )
    colsList[["cds"]] <- c(colsList[["cds"]], colsList[["gene"]])
    colsList[["exon"]] <- c(colsList[["exon"]], colsList[["gene"]])
    colsList[["tx"]] <- c(colsList[["tx"]], colsList[["gene"]])
    colsList <- lapply(
        X = colsList,
        FUN = function(x) {
            x <- sort(unique(tolower(x)))
            gsub(
                pattern = "^(cds|exon|gene|tx)",
                replacement = "\\1_",
                x = x
            )
        }
    )
    args <- list(
        "x" = object,
        "columns" = colsList[[level]]
    )
    switch(
        EXPR = level,
        "genes" = {
            args <- append(
                x = args,
                values = list(
                    "single.strand.genes.only" = TRUE
                )
            )
        }
    )
    what <- get(
        x = level,
        envir = asNamespace("GenomicFeatures"),
        inherits = FALSE
    )
    assert(is.function(what))
    suppressMessages({
        gr <- do.call(what = what, args = args)
    })
    ## This will also return metadata slotted into `genomeInfo`.
    meta <- metadata(gr)
    gffMeta <- attr(x = txdb, which = "gffMetadata", exact = TRUE)
    if (is.list(gffMeta)) {
        meta <- append(x = meta, values = gffMeta)
    }
    meta[["level"]] <- level
    metadata(gr) <- meta
    gr <- .makeGRanges(
        object = gr,
        ignoreVersion = ignoreVersion
    )
    gr
}
