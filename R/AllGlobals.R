#' All global variables
#' @noRd
NULL



#' GFF file name pattern matching
#'
#' @note Updated 2021-01-23.
#' @noRd
.gffPatterns <- list(
    "ensembl" = paste0(
        ## e.g. "Homo_sapiens.GRCh38.102.gtf.gz"
        "^([a-z0-9]+_)?",              # BiocFileCache.
        "([A-Z][a-z]+_[a-z]+)",        # "Homo_sapiens".
        "\\.([A-Za-z0-9]+)",           # "GRCh38"
        "\\.([0-9]+)",                 # "102".
        "(\\.chr_patch_hapl_scaff)?",
        "\\.(gff3|gtf)",
        "(\\.gz)?$"
    ),
    "flybase" = paste0(
        ## e.g. "dmel-all-r6.37.gtf.gz".
        "^([a-z0-9]+_)?",       # BiocFileCache
        "^([^-]+)",             # "dmel"
        "-([^-]+)",             # "all"
        "-(r[0-9]+\\.[0-9]+)",  # "r6.37"
        "\\.(gff|gtf)",
        "(\\.gz)?$"
    ),
    "gencode" = paste0(
        ## e.g. "gencode.v36.annotation.gtf.gz"
        "^([a-z0-9]+_)?",  # BiocFileCache.
        "gencode",
        "\\.v([M0-9]+)",   # "36" (human) / "M25" (mouse)
        "(lift37)?",       # GRCh37-specific.
        "\\.annotation",
        "\\.(gff3|gtf)",
        "(\\.gz)?$"
    ),
    "ucsc" = paste0(
        ## e.g. "hg38.ensGene.gtf.gz"
        "^([0-9a-z]_)?",                           # BiocFileCache.
        "([a-z]+[A-Za-z]+[0-9]+)",                 # "hg38".
        "(ensGene|knownGene|ncbiRefSeq|refGene)",  # "ensGene"
        "\\.gtf",
        "(\\.gz)?$"
    ),
    "wormbase" = paste0(
        ## e.g. c_elegans.PRJNA13758.WS279.canonical_geneset.gtf.gz
        "^([a-z0-9]+_)?",   # BiocFileCache.
        "^([a-z]_[a-z]+)",  # "c_elegans".
        "\\.([A-Z0-9]+)",   # "PRJNA13758".
        "\\.(WS[0-9]+)",    # "WS279".
        "\\.([a-z_]+)",     # "canonical_geneset".
        "\\.(gff3|gtf)",
        "(\\.gz)?$"
    )
)



#' Package version
#'
#' @note Updated 2020-10-06.
#' @noRd
.version <- packageVersion(packageName())



#' AcidGenomes test data URL
#'
#' @export
#' @keywords internal
#' @note Updated 2020-10-06.
#'
#' @examples
#' AcidGenomesTestsURL
AcidGenomesTestsURL <-  # nolint
    paste0(
        "https://tests.acidgenomics.com/AcidGenomes/",
        "v", .version$major, ".", .version$minor  # nolint
    )
