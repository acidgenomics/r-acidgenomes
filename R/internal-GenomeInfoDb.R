#' Get Seqinfo
#'
#' @note Updated 2022-02-08.
#' @noRd
#'
#' @param x GFF file or `getGFFMetadata()` return list.
#'
#' @return `Seqinfo` or `NULL` (on failure).
#'
#' @examples
#' ## Ensembl ====
#' file <- pasteURL(
#'     "ftp.ensembl.org",
#'     "pub",
#'     "release-102",
#'     "gtf",
#'     "homo_sapiens",
#'     "Homo_sapiens.GRCh38.102.gtf.gz",
#'     protocol = "ftp"
#' )
#' seq <- .getSeqinfo(file)
#' print(seq)
#'
#' ## GENCODE ====
#' file <- pasteURL(
#'     "ftp.ebi.ac.uk",
#'     "pub",
#'     "databases",
#'     "gencode",
#'     "Gencode_human",
#'     "release_36",
#'     "gencode.v36.annotation.gtf.gz",
#'     protocol = "ftp"
#' )
#' seq <- .getSeqinfo(file)
#' print(seq)
#'
#' ## RefSeq ====
#' file <- pasteURL(
#'     "ftp.ncbi.nlm.nih.gov",
#'     "genomes",
#'     "refseq",
#'     "vertebrate_mammalian",
#'     "Homo_sapiens",
#'     "all_assembly_versions",
#'     "GCF_000001405.38_GRCh38.p12",
#'     "GCF_000001405.38_GRCh38.p12_genomic.gff.gz",
#'     protocol = "ftp"
#' )
#' seq <- .getSeqinfo(file)
#' print(seq)
#'
#' ## UCSC ====
#' file <- pasteURL(
#'     "hgdownload.soe.ucsc.edu",
#'     "goldenPath",
#'     "hg38",
#'     "bigZips",
#'     "genes",
#'     "hg38.ensGene.gtf.gz",
#'     protocol = "ftp"
#' )
#' seq <- .getSeqinfo(file)
#' print(seq)
.getSeqinfo <- function(x) {
    pkgs <- .packages()
    ## Allowing pass-in of either file or GFF metadata list here.
    if (!is.list(x)) {
        x <- .getGFFMetadata(x)
    }
    assert(is.list(x))
    if (
        !isSubset(
            x = c("genomeBuild", "organism", "provider"),
            y = names(x)
        ) ||
        any(
            is.na(c(
                x[["genomeBuild"]],
                x[["organism"]],
                x[["provider"]]
            ))
        )
    ) {
        return(NULL)
    }
    assert(
        isString(x[["file"]]),
        isString(x[["genomeBuild"]]),
        isString(x[["organism"]]),
        isString(x[["provider"]]),
        isScalar(x[["release"]]) || is.null(x[["release"]]),
        msg = "GFF file does not contain sufficient metadata."
    )
    seq <- NULL
    tryCatch(
        expr = suppressPackageStartupMessages({
            switch(
                EXPR = x[["provider"]],
                "Ensembl" = {
                    assert(isInt(x[["release"]]))
                    args <- list(
                        "species" = x[["organism"]],
                        "release" = x[["release"]],
                        "as.Seqinfo" = TRUE
                    )
                    ## The `use.grch37` flag isn't currently working with
                    ## GenomeInfoDb v1.26.2, but may be improved in the future.
                    if (isMatchingFixed(
                        pattern = "GRCh37",
                        x = x[["genomeBuild"]])
                    ) {
                        args[["use.grch37"]] <- TRUE
                    }
                    seq <- do.call(
                        what = getChromInfoFromEnsembl,
                        args = args
                    )
                },
                "GENCODE" = {
                    genome <- .mapGenomeBuildToUCSC(x[["genomeBuild"]])
                    seq <- Seqinfo(genome = genome)
                    genome(seq) <- x[["genomeBuild"]]
                },
                "RefSeq" = {
                    seq <- .getRefSeqSeqinfo(file)
                },
                "UCSC" = {
                    seq <- Seqinfo(genome = x[["genomeBuild"]])
                }
            )
        }),
        error = function(e) {
            alertWarning(sprintf(
                "Automatic {.cls %s} assignment failed.", "Seqinfo"
            ))
            NULL
        }
    )
    assert(isAny(seq, c("Seqinfo", "NULL")))
    if (is(seq, "Seqinfo")) {
        assert(
            identical(
                x = unique(unname(genome(seq))),
                y = x[["genomeBuild"]]
            ),
            msg = sprintf(
                "Seqinfo does not contain expected genome build: '%s'.",
                x[["genomeBuild"]]
            )
        )
    }
    forceDetach(keep = pkgs)
    seq
}
