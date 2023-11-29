#' Get Seqinfo
#'
#' @note Updated 2023-11-29.
#' @noRd
#'
#' @param x GFF file or `.getGffMetadata()` return list.
#'
#' @return `Seqinfo` or `NULL` (on failure).
#'
#' @examples
#' ## Ensembl ====
#' file <- pasteUrl(
#'     "ftp.ensembl.org",
#'     "pub",
#'     "release-102",
#'     "gtf",
#'     "homo_sapiens",
#'     "Homo_sapiens.GRCh38.102.gtf.gz",
#'     protocol = "https"
#' )
#' seq <- .getSeqinfo(file)
#' print(seq)
#'
#' ## GENCODE ====
#' file <- pasteUrl(
#'     "ftp.ebi.ac.uk",
#'     "pub",
#'     "databases",
#'     "gencode",
#'     "Gencode_human",
#'     "release_36",
#'     "gencode.v36.annotation.gtf.gz",
#'     protocol = "https"
#' )
#' seq <- .getSeqinfo(file)
#' print(seq)
#'
#' ## RefSeq ====
#' file <- pasteUrl(
#'     "ftp.ncbi.nlm.nih.gov",
#'     "genomes",
#'     "refseq",
#'     "vertebrate_mammalian",
#'     "Homo_sapiens",
#'     "all_assembly_versions",
#'     "GCF_000001405.38_GRCh38.p12",
#'     "GCF_000001405.38_GRCh38.p12_genomic.gff.gz",
#'     protocol = "https"
#' )
#' seq <- .getSeqinfo(file)
#' print(seq)
#'
#' ## UCSC ====
#' file <- pasteUrl(
#'     "hgdownload.soe.ucsc.edu",
#'     "goldenPath",
#'     "hg38",
#'     "bigZips",
#'     "genes",
#'     "hg38.knownGene.gtf.gz",
#'     protocol = "ftp"
#' )
#' seq <- .getSeqinfo(file)
#' print(seq)
.getSeqinfo <- function(x) {
    if (!is.list(x)) {
        x <- .getGffMetadata(x)
    }
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
    seq <- switch(
        EXPR = x[["provider"]],
        "Ensembl" = {
            .getEnsemblSeqinfo(
                organism = x[["organism"]],
                genomeBuild = x[["genomeBuild"]],
                release = x[["release"]]
            )
        },
        "GENCODE" = {
            .getGencodeSeqinfo(
                organism = x[["organism"]],
                genomeBuild = x[["genomeBuild"]],
                release = x[["release"]]
            )
        },
        "RefSeq" = {
            .getRefSeqSeqinfo(
                file = ifelse(
                    test = isAUrl(x[["url"]]),
                    yes = x[["url"]],
                    no = x[["file"]]
                )
            )
        },
        "UCSC" = {
            .getUcscSeqinfo(genomeBuild = x[["genomeBuild"]])
        },
        NULL
    )
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
        seq <- seq[sort(seqnames(seq))]
    }
    seq
}



#' Get Ensembl genome assembly seqinfo
#'
#' @note Updated 2023-11-22.
#' @noRd
.getEnsemblSeqinfo <- function(organism, genomeBuild, release) {
    assert(
        isString(organism),
        isString(genomeBuild),
        isInt(release)
    )
    args <- list("species" = organism, "release" = release, "as.Seqinfo" = TRUE)
    if (grepl(pattern = "GRCh37", x = genomeBuild, fixed = TRUE)) {
        args[["use.grch37"]] <- TRUE
    }
    quietly({
        seq <- do.call(what = getChromInfoFromEnsembl, args = args)
    })
    assert(is(seq, "Seqinfo"))
    seq
}



#' Get GENCODE genome assembly seqinfo
#'
#' @note Updated 2023-11-22.
#' @noRd
.getGencodeSeqinfo <- function(organism, genomeBuild, release) {
    seq <- .getEnsemblSeqinfo(
        organism = organism,
        genomeBuild = genomeBuild,
        release = mapGencodeToEnsembl(release)
    )
    assert(is(seq, "Seqinfo"))
    genome(seq) <- sub(
        pattern = "\\.p[0-9]+$",
        replacement = "",
        x = genome(seq)
    )
    keep <- intersect(
        x = seqnames(seq),
        y = c(seq(from = 1L, to = 23L), "MT", "X", "Y")
    )
    assert(hasLength(keep))
    seq <- seq[keep]
    seqnames(seq)[seqnames(seq) == "MT"] <- "M"
    seqnames(seq) <- paste0("chr", seqnames(seq))
    seq
}



#' Get UCSC genome assembly seqinfo
#'
#' @note Updated 2023-01-30.
#' @noRd
.getUcscSeqinfo <- function(genomeBuild) {
    tryCatch(
        expr = {
            quietly({
                Seqinfo(genome = genomeBuild)
            })
        },
        error = function(e) {
            NULL
        }
    )
}
