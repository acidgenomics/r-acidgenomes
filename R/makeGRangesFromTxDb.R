#' Make GRanges from TxDb object
#'
#' @export
#' @note Updated 2021-01-26.
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
    ignoreVersion = FALSE,
    synonyms = FALSE
) {
    assert(is(object, "TxDb"))
    level <- match.arg(level)
    keys <- columns(object)
    colsList <- list(
        "cds" = grep(pattern = "^CDS", x = keys, value = TRUE),
        "exons" = grep(pattern = "^EXON", x = keys, value = TRUE),
        "genes" = grep(pattern = "^GENE", x = keys, value = TRUE),
        "transcripts" = grep(pattern = "^TX", x = keys, value = TRUE)
    )
    colsList[["cds"]] <-
        c(colsList[["cds"]], colsList[["genes"]])
    colsList[["exons"]] <-
        c(colsList[["exons"]], colsList[["genes"]])
    colsList[["transcripts"]] <-
        c(colsList[["transcripts"]], colsList[["genes"]])
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
    columns <- colsList[[level]]
    assert(isCharacter(columns))
    args <- list(
        "x" = object,
        "columns" = columns
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
    if (isSubset(c("tx_id", "tx_name"), colnames(mcols(gr)))) {
        ## Improve identifier handling for RefSeq input.
        if (is.integer(decode(mcols(gr)[["tx_id"]]))) {
            mcols(gr)[["tx_number"]] <- mcols(gr)[["tx_id"]]
            mcols(gr)[["tx_id"]] <- mcols(gr)[["tx_name"]]
        }
    }
    ## This will also return metadata slotted into `genomeInfo`.
    meta <- metadata(gr)
    gffMeta <- attr(x = object, which = "gffMetadata", exact = TRUE)
    if (is.list(gffMeta)) {
        meta <- append(x = meta, values = gffMeta)
    }
    meta[["level"]] <- level
    metadata(gr) <- meta
    gr <- .makeGRanges(
        object = gr,
        ignoreVersion = ignoreVersion,
        synonyms = synonyms
    )
    gr
}
