#' Make genomic ranges (`GRanges`) from TxDb object
#'
#' @note Updated 2023-04-26.
#' @noRd
#'
#' @inheritParams AcidRoxygen::params
#' @inheritParams params
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
#' txdb <- .makeTxDbFromGFF(file)
#' gr <- .makeGRangesFromTxDb(object = txdb, ignoreVersion = FALSE)
#' print(gr)
.makeGRangesFromTxDb <-
    function(object,
             level = c("transcripts", "genes", "exons", "cds"),
             ignoreVersion = TRUE) {
        assert(
            requireNamespaces("AnnotationDbi"),
            is(object, "TxDb")
        )
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
            ## Ensure we coerce gene identifiers to character vector. This
            ## currently gets returned as CharacterList for RefSeq.
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
            ## Improve identifier handling for UCSC and/or RefSeq input. Note
            ## that RefSeq transcript names currently map to the gene names
            ## here, which is incorrect and confusing.
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
            ## Drop any remaining elements where the transcript and gene
            ## identifiers are identical. This is garbage output from RefSeq
            ## that we don't want to include in transcript-to-gene mappings.
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
            ignoreVersion = ignoreVersion
        )
        gr
    }



## nolint start

#' Make TxDb from a GFF/GTF file
#'
#' Wrapper for GenomicFeatures `makeTxDbFromGFF` importer.
#'
#' @note Updated 2022-01-12.
#' @noRd
#'
#' @return `TxDb`.
#'
#' @seealso
#' - `GenomicFeatures::makeTxDbFromGFF()`.
#' - `GenomicFeatures::supportedMiRBaseBuildValues()`.
#' Note that *Homo sapiens* GRCh38 isn't currently supported in mirbase.db.
#' - `TxDb.Hsapiens.UCSC.hg38.knownGene` package.`
#'
#' @examples
#' ## GENCODE ====
#' ## > gtfFile <- pasteURL(
#' ## >     "ftp.ebi.ac.uk",
#' ## >     "pub",
#' ## >     "databases",
#' ## >     "gencode",
#' ## >     "Gencode_human",
#' ## >     "release_36",
#' ## >     "gencode.v36.annotation.gtf.gz",
#' ## >     protocol = "ftp"
#' ## > )
#' ## > txdb <- .makeTxDbFromGFF(gtfFile)
#' ## > print(txdb)
#'
#' ## RefSeq ====
#' ## > gffFile <- pasteURL(
#' ## >     "ftp.ncbi.nlm.nih.gov",
#' ## >     "genomes",
#' ## >     "refseq",
#' ## >     "vertebrate_mammalian",
#' ## >     "Homo_sapiens",
#' ## >     "all_assembly_versions",
#' ## >     "GCF_000001405.38_GRCh38.p12",
#' ## >     "GCF_000001405.38_GRCh38.p12_genomic.gff.gz",
#' ## >     protocol = "ftp"
#' ## > )
#' ## > txdb <- .makeTxDbFromGFF(gffFile)
#' ## > print(txdb)
NULL

## nolint end

.makeTxDbFromGFF <- function(file, meta) {
    assert(
        .isSupportedGFF(file),
        is.list(meta)
    )
    ## Check for input of unsupported files.
    ## See `.gffPatterns` for details.
    denylist <- c(
        "refseq_gtf" = paste0(
            "^([0-9a-z]_)?",
            "(GC[AF]_[0-9]+\\.[0-9]+)",
            "_([^_]+)",
            "_(.+)",
            "\\.gtf",
            "(\\.gz)?$"
        )
    )
    if (grepl(
        pattern = denylist[["refseq_gtf"]],
        x = basename(file)
    )) {
        abort(sprintf(
            paste(
                "Unsupported file: {.file %s}.",
                "Use RefSeq GFF instead of GTF.",
                sep = "\n"
            ),
            basename(file)
        ))
    }
    alert(sprintf(
        "Making {.cls %s} from {.file %s} with {.pkg %s}::{.fun %s}.",
        "TxDb", file,
        "GenomicFeatures", "makeTxDbFromGFF"
    ))
    assert(requireNamespaces("GenomicFeatures"))
    if (isAFile(file)) {
        file <- realpath(file)
    }
    seqinfo <- .getSeqinfo(meta)
    assert(isAny(seqinfo, c("Seqinfo", "NULL")))
    ## Additional arguments of potential future interest:
    ## - dbxrefTag: This can help override primary identifier to use.
    ## - miRBaseBuild: miRBase annotations could be useful for future genome
    ## builds. Note that it's currently out of date with GRCh38.
    ## https://github.com/Bioconductor/GenomicFeatures/issues/27
    args <- list(
        "file" = .cacheIt(file),
        "dataSource" = file,
        "organism" = meta[["organism"]]
    )
    if (!is.null(seqinfo)) {
        args <- append(x = args, values = list("chrominfo" = seqinfo))
    }
    what <- GenomicFeatures::makeTxDbFromGFF
    suppressWarnings({
        txdb <- do.call(what = what, args = args)
    })
    assert(is(txdb, "TxDb"))
    ## Stash the GFF metadata, so we can access in `makeGRangesFromGFF()`.
    attr(txdb, which = "gffMetadata") <- meta
    assert(validObject(txdb))
    txdb
}
