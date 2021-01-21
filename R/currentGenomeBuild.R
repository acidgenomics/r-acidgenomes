#' Current genome build
#'
#' Fetch the current genome build (assembly) version from online resources.
#'
#' @name currentGenomeBuild
#' @note Updated 2020-01-14.
#'
#' @inheritParams AcidRoxygen::params
#'
#' @return `character(1)`.
#'   Genome assembly build version.
#'
#' @seealso
#' - [Ensembl REST API](https://rest.ensembl.org).
#' - [UCSC REST API](https://genome.ucsc.edu/goldenPath/help/api.html).
#' - [RefSeq genomes](https://ftp.ncbi.nlm.nih.gov/genomes/refseq/).
#'
#' @examples
#' currentEnsemblBuild("Homo sapiens")
#' currentGencodeBuild("Homo sapiens")
#' currentRefSeqGenomeBuild(
#'     organism = "Homo sapiens",
#'     taxonomicGroup = "vertebrate_mammalian"
#' )
#' currentUCSCGenomeBuild("Homo sapiens")
NULL



#' Return a simple (minimal) genome build version
#'
#' @details
#' For example, sanitize "GRCh38.p13" to simply "GRCh38".
#'
#' @note Updated 2021-01-07.
#' @noRd
.simpleGenomeBuild <- function(x) {
    assert(isString(x))
    x <- sub(pattern = "\\.[^\\.]+$", replacement = "", x = x)
    x
}



## Updated 2021-01-21.
#' @rdname currentGenomeBuild
#' @export
currentEnsemblGenomeBuild <- function(organism) {
    assert(isOrganism(organism))
    organism <- snakeCase(organism)
    json <- getJSON(pasteURL(
        "rest.ensembl.org",
        "info",
        "assembly",
        paste0(organism, "?"),
        protocol = "https"
    ))
    assert(isSubset("assembly_name", names(json)))
    out <- json[["assembly_name"]]
    assert(isString(out))
    out
}



## Updated 2021-01-05.
#' @rdname currentGenomeBuild
#' @export
currentGencodeGenomeBuild <- function(organism) {
    organism <- match.arg(
        arg = organism,
        choices = c("Homo sapiens", "Mus musculus")
    )
    currentEnsemblBuild(organism)
}



## Alternate approach using URL only:
## https://ftp.ncbi.nlm.nih.gov/genomes/refseq/<taxonomic_group>/<organism>/
##     latest_assembly_versions/

## Updated 2021-01-14.
#' @rdname currentGenomeBuild
#' @param taxonomicGroup `character(1)`.
#'   *Only applies to RefSeq*.
#'   FTP server taxonomic group subdirectory path (e.g. "vertebrate_mammalian").
#'   Defining this manually avoids having to query the FTP server.
#' @export
currentRefSeqGenomeBuild <- function(
    organism,
    taxonomicGroup = NULL
) {
    assert(
        isOrganism(organism),
        isString(taxonomicGroup, nullOK = TRUE)
    )
    baseURL <- .getRefSeqGenomeURL(
        organism = organism,
        taxonomicGroup = taxonomicGroup,
        quiet = TRUE
    )
    summary <- getRefSeqAssemblySummary(
        file = pasteURL(baseURL, "assembly_summary.txt")
    )
    assert(isSubset("ftp_path", names(summary)))
    out <- basename(summary[["ftp_path"]])
    out
}



## Updated 2021-01-05.
#' @rdname currentGenomeBuild
#' @export
currentUCSCGenomeBuild <- function(organism) {
    assert(isOrganism(organism))
    json <- getJSON("https://api.genome.ucsc.edu/list/ucscGenomes")
    assert(isSubset("ucscGenomes", names(json)))
    json <- json[["ucscGenomes"]]
    l <- mapply(
        name = names(json),
        x = json,
        FUN = function(name, x) {
            ## Other useful keys: description, sourceName.
            c(
                "build" = name,
                "active" = x[["active"]],
                "orderKey" = x[["orderKey"]],
                "scientificName" = x[["scientificName"]]
            )
        },
        SIMPLIFY = FALSE,
        USE.NAMES = FALSE
    )
    df <- map_df(.x = l, .f = unlist)
    df <- df[df[["active"]] == 1L, , drop = FALSE]
    if (!isSubset(organism, unique(df[["scientificName"]]))) {
        stop("Invalid organism.")
    }
    df <- df[df[["scientificName"]] == organism, , drop = FALSE]
    ## The latest genome build has the lowest "orderKey" value.
    df <- df[order(df[["orderKey"]]), , drop = FALSE]
    out <- df[["build"]][[1L]]
    assert(isString(out))
    out
}
