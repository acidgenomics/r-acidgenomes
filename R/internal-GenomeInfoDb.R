## FIXME Getting the Seqinfo from UCSC now appears to be broken on Bioconductor.
## Error in .order_seqlevels(chrom_sizes[, "chrom"]) :
##   !anyNA(m32) is not TRUE
## Calls: <Anonymous> ... FETCH_ORDERED_CHROM_SIZES -> .order_seqlevels -> stopifnot
## Backtrace:
##     ▆
##  1. └─GenomeInfoDb::Seqinfo(genome = "hg38")
##  2.   └─GenomeInfoDb:::.make_Seqinfo_from_genome(genome)
##  3.     └─GenomeInfoDb::getChromInfoFromUCSC(genome, as.Seqinfo = TRUE)
##  4.       └─GenomeInfoDb:::.get_chrom_info_for_registered_UCSC_genome(...)
##  5.         └─GenomeInfoDb:::.get_raw_chrom_info_for_registered_UCSC_genome(...)
##  6.           └─GenomeInfoDb:::.fetch_raw_chrom_info_from_UCSC(...)
##  7.             └─GenomeInfoDb (local) FETCH_ORDERED_CHROM_SIZES(goldenPath.url = goldenPath.url)
##  8.               └─GenomeInfoDb (local) .order_seqlevels(chrom_sizes[, "chrom"]) at registered/UCSC_genomes/hg38.R:69:4
##  9.                 └─base::stopifnot(!anyNA(m32)) at registered/UCSC_genomes/hg38.R:37:4



#' Get Seqinfo
#'
#' @note Updated 2023-01-30.
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
            anyNA(c(
                x[["genomeBuild"]],
                x[["organism"]],
                x[["provider"]]
            ))
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
                    seq <- .getEnsemblSeqinfo(
                        organism = x[["organism"]],
                        genomeBuild = x[["genomeBuild"]],
                        release = x[["release"]]
                    )
                },
                "GENCODE" = {
                    seq <- .getEnsemblSeqinfo(
                        organism = x[["organism"]],
                        genomeBuild = x[["genomeBuild"]],
                        release = mapGencodeToEnsembl(x[["release"]])
                    )
                },
                "RefSeq" = {
                    seq <- .getRefSeqSeqinfo(
                        file = ifelse(
                            test = isAURL(x[["url"]]),
                            yes = x[["url"]],
                            no = x[["file"]]
                        )
                    )
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



## Updated 2023-01-30.
.getEnsemblSeqinfo <- function(organism, genomeBuild, release) {
    assert(
        isString(organism),
        isString(genomeBuild),
        isInt(release)
    )
    args <- list(
        "species" = organism,
        "release" = release,
        "as.Seqinfo" = TRUE
    )
    ## The `use.grch37` flag isn't currently working with
    ## GenomeInfoDb v1.26.2, but may be improved in the future.
    if (grepl(
        pattern = "GRCh37",
        x = genomeBuild,
        fixed = TRUE
    )
    ) {
        args[["use.grch37"]] <- TRUE
    }
    do.call(
        what = getChromInfoFromEnsembl,
        args = args
    )
}
