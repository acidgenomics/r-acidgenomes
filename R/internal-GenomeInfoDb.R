#' Get Seqinfo
#'
#' @note Updated 2021-02-01.
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
    if (!is.list(x)) {
        x <- getGFFMetadata(x)
    }
    assert(is.list(x))
    if (!isSubset(x = c("genomeBuild", "organism", "provider"), y = names(x))) {
        return(NULL)
    }
    assert(
        isString(x[["file"]]),
        isString(x[["genomeBuild"]]),
        isString(x[["organism"]]),
        isString(x[["provider"]]),
        isScalar(x[["release"]]) || is.null(x[["release"]])
    )
    seq <- tryCatch(
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
                    ## NOTE This step is failing for Ensembl 103 with
                    ## GenomeInfoDb 1.26.2.
                    do.call(what = getChromInfoFromEnsembl, args = args)
                },
                "GENCODE" = {
                    Seqinfo(genome = mapNCBIBuildToUCSC(x[["genomeBuild"]]))
                },
                "RefSeq" = {
                    .getRefSeqSeqinfo(x[["file"]])
                },
                "UCSC" = {
                    Seqinfo(genome = x[["genomeBuild"]])
                }
            )
        }),
        error = function(e) {
            alertWarning("Automatic {.var seqinfo} assignment failed.")
            NULL
        }
    )
    assert(isAny(seq, c("Seqinfo", "NULL")))
    forceDetach(keep = pkgs)
    seq
}
