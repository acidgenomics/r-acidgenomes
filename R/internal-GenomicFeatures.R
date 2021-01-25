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
    level = c("transcripts", "genes", "exons", "cds")
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
    ## This will return metadata slotted into `genomeInfo`.
    ##  [1] "Db type"
    ##  [2] "Supporting package"
    ##  [3] "Data source"
    ##  [4] "Organism"
    ##  [5] "Taxonomy ID"
    ##  [6] "miRBase build ID"
    ##  [7] "Genome"
    ##  [8] "Nb of transcripts"
    ##  [9] "Db created by"
    ## [10] "Creation time"
    ## [11] "GenomicFeatures version at creation time"
    ## [12] "RSQLite version at creation time"
    ## [13] "DBSCHEMAVERSION"

    idCol <- .matchGRangesNamesColumn(gr)
    assert(is(gr, "GRanges"))
    gr
}
