#' Get Seqinfo
#'
#' @note Updated 2021-01-25.
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
    if (!is.list(x)) {
        x <- getGFFMetadata(x)
    }
    assert(
        is.list(x),
        isSubset(
            x = c("genomeBuild", "organism", "provider"),
            y = names(x)
        )
    )
    genomeBuild <- x[["genomeBuild"]]
    organism <- x[["organism"]]
    provider <- x[["provider"]]
    release <- x[["release"]]
    assert(
        isString(genomeBuild),
        isString(organism),
        isString(provider),
        isScalar(release) || is.null(release)
    )
    out <- tryCatch(
        expr = {
            genome <- switch(
                EXPR = provider,
                "Ensembl" = {
                    ## This may support a `use.grch37` flag, but it's not
                    ## currently tested well according to documentation.
                    if (isMatchingFixed(pattern = "GRCh37", x = genomeBuild)) {
                        return(NULL)
                    }
                    getChromInfoFromEnsembl(
                        species = organism,
                        release = release,
                        as.Seqinfo = TRUE
                    )
                },
                "GENCODE" = {
                    genome <- unname(mapNCBIBuildToUCSC(genomeBuild))
                    Seqinfo(genome = genome)
                },
                "RefSeq" = {
                    .getRefSeqSeqinfo(file)
                },
                "UCSC" = {
                    genome <- genomeBuild
                    Seqinfo(genome = genome)
                },
                NULL
            )

        },
        error = function(e) {
            alertWarning(paste(
                "Automatic seqinfo detection failed.",
                "Manual input of {.var seqinfo} is recommended."
            ))
            NULL
        }
    )
    assert(isAny(out, c("Seqinfo", "NULL")))
    out
}
