#' Make GenomicRanges from TxDb object
#'
#' @export
#' @note Updated 2021-04-27.
#'
#' @inheritParams AcidRoxygen::params
#' @inheritParams params
#'
#' @return `GenomicRanges`.
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
#' txdb <- AcidGenomes::makeTxDbFromGFF(file)
#' gr <- makeGRangesFromTxDb(object = txdb, ignoreVersion = FALSE)
#' print(gr)
makeGRangesFromTxDb <- function(
    object,
    level = c("transcripts", "genes", "exons", "cds"),
    ignoreVersion = TRUE,
    synonyms = FALSE
) {
    pkgs <- .packages()
    requireNamespaces("AnnotationDbi")
    assert(is(object, "TxDb"))
    level <- match.arg(level)
    cols <- AnnotationDbi::columns(object)
    colsList <- list(
        "cds" = grep(pattern = "^CDS", x = cols, value = TRUE),
        "exons" = grep(pattern = "^EXON", x = cols, value = TRUE),
        "genes" = grep(pattern = "^GENE", x = cols, value = TRUE),
        "transcripts" = grep(pattern = "^TX", x = cols, value = TRUE)
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
    ## Transcript-specific fixes.
    if (identical(level, "transcripts")) {
        ## Ensure we coerce gene identifiers to character vector, if necessary.
        ## This currently gets returned as CharacterList for RefSeq.
        if (is(mcols(gr)[["gene_id"]], "CharacterList")) {
            mcols(gr)[["gene_id"]] <-
                as.character(mcols(gr)[["gene_id"]])
        }
        ## Drop transcripts that don't map to genes.
        keep <- !is.na(mcols(gr)[["gene_id"]])
        gr <- gr[keep]
        ## Ensure "rna-" prefix is correctly removed from identifiers.
        ## This is not currently handled correctly for RefSeq input.
        ## (e.g. "rna-MIR1302-2", "rna-TRNP", etc.).
        if (any(grepl(pattern = "^rna-", x = mcols(gr)[["tx_id"]]))) {
            mcols(gr)[["tx_id"]] <-
                gsub(
                    pattern = "^rna-",
                    replacement = "",
                    x = mcols(gr)[["tx_id"]]
                )
        }
        if (any(grepl(pattern = "^rna-", x = mcols(gr)[["tx_name"]]))) {
            mcols(gr)[["tx_name"]] <-
                gsub(
                    pattern = "^rna-",
                    replacement = "",
                    x = mcols(gr)[["tx_name"]]
                )
        }
        ## Improve identifier handling for UCSC and/or RefSeq input. Note that
        ## RefSeq transcript names currently map to the gene names here, which
        ## is incorrect and confusing.
        if (
            isSubset(c("tx_id", "tx_name"), colnames(mcols(gr))) &&
            is.integer(decode(mcols(gr)[["tx_id"]]))
        ) {
            ## Not sure these numbers are actually useful, but keep for the
            ## time being just in case.
            mcols(gr)[["tx_number"]] <- mcols(gr)[["tx_id"]]
            mcols(gr)[["tx_id"]] <- mcols(gr)[["tx_name"]]
        }

        ## Drop any transcript identifiers that return NA. This can happen
        ## with RefSeq return.
        keep <- !is.na(mcols(gr)[["tx_id"]])
        gr <- gr[keep]
        ## Drop any remaining elements where the transcript and gene identifiers
        ## are identical. This is garbage output from RefSeq that we don't want
        ## to include in transcript-to-gene mappings.
        keep <- apply(
            X = mcols(gr),
            MARGIN = 1L,
            FUN = function(x) {
                !identical(
                    x = as.character(x[["tx_id"]]),
                    x[["gene_id"]]
                )
            }
        )
        gr <- gr[keep]
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
    forceDetach(keep = pkgs)
    gr
}
