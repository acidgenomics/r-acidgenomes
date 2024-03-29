#' All global variables
#' @noRd
NULL



#' GFF file name pattern matching
#'
#' @note Updated 2021-02-08.
#' @noRd
.gffPatterns <- list(
    "ensembl" = paste0(
        ## "Homo_sapiens.GRCh38.102.gtf.gz"
        "^([a-z0-9]+_)?", # BiocFileCache
        "([A-Z][a-z]+_[a-z]+)", # "Homo_sapiens"
        "\\.([A-Za-z0-9]+)", # "GRCh38"
        "\\.([0-9]+)", # "102"
        "(\\.chr_patch_hapl_scaff)?",
        "\\.(gff3|gtf)",
        "(\\.gz)?$"
    ),
    "flybase" = paste0(
        ## "dmel-all-r6.37.gtf.gz".
        "^([a-z0-9]+_)?", # BiocFileCache
        "^([^-]+)", # "dmel"
        "-([^-]+)", # "all"
        "-(r[0-9]+\\.[0-9]+)", # "r6.37"
        "\\.(gff|gtf)",
        "(\\.gz)?$"
    ),
    "gencode" = paste0(
        ## "gencode.v36.annotation.gtf.gz"
        "^([a-z0-9]+_)?", # BiocFileCache
        "gencode",
        "\\.v([M0-9]+)", # "36" (human); "M25" (mouse)
        "(lift37)?", # GRCh37-specific
        "\\.annotation",
        "\\.(gff3|gtf)",
        "(\\.gz)?$"
    ),
    "refseq" = paste0(
        ## "GCF_000001405.38_GRCh38.p12_genomic.gff.gz"
        ## "GCA_000001405.15_GRCh38_full_analysis_set.refseq_annotation.gff.gz"
        "^([0-9a-z]+_)?", # BiocFileCache
        "(GC[AF]_[0-9]+\\.[0-9]+)", # "GCF_000001405.38"
        "_([^_]+)", # "GRCh38.p12"
        "_(.+)", # "genomic" or "full_analysis_set"
        "\\.(gff|gtf)",
        "(\\.gz)?$"
    ),
    "ucsc" = paste0(
        ## "hg38.knownGene.gtf.gz"
        "^([0-9a-z]+_)?", # BiocFileCache
        "([a-z]+[A-Za-z]+[0-9]+)", # "hg38"
        "\\.(ensGene|knownGene|ncbiRefSeq|refGene)", # "ensGene"
        "\\.gtf",
        "(\\.gz)?$"
    ),
    "wormbase" = paste0(
        ## "c_elegans.PRJNA13758.WS279.canonical_geneset.gtf.gz"
        "^([a-z0-9]+_)?", # BiocFileCache
        "^([a-z]_[a-z]+)", # "c_elegans"
        "\\.([A-Z0-9]+)", # "PRJNA13758"
        "\\.(WS[0-9]+)", # "WS279"
        "\\.([a-z_]+)", # "canonical_geneset"
        "\\.(gff3|gtf)",
        "(\\.gz)?$"
    )
)



#' GRanges annotation levels
#'
#' @note Updated 2023-04-26.
#' @noRd
.grangesLevels <- c(
    "cds",
    "exons",
    "genes",
    "transcripts"
)



#' Package name
#'
#' @note Updated 2021-03-03.
#' @noRd
.pkgName <- packageName()



#' Package version
#'
#' @note Updated 2021-03-03.
#' @noRd
.pkgVersion <- packageVersion(.pkgName)



#' AcidGenomes test data URL
#'
#' @export
#' @keywords internal
#' @note Updated 2023-10-04.
#'
#' @examples
#' AcidGenomesTestsUrl
AcidGenomesTestsUrl <- # nolint
    "https://r.acidgenomics.com/testdata/acidgenomes"
