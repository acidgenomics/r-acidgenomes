#' Current genome build
#'
#' Fetch the current genome build (assembly) version from online resources.
#'
#' @name currentGenomeBuild
#' @note Updated 2020-01-05.
#'
#' @inheritParams AcidRoxygen::params
#' @param subdir `character(1)`.
#'   *Only applies to RefSeq*.
#'   FTP server subdirectory path (e.g. "vertebrate_mammalian").
#'
#' @return `character(1)`.
#'   Genome assembly build version.
#'
#' @seealso
#' - [Ensembl REST API](https://rest.ensembl.org).
#' - [UCSC Genome Browser REST API](https://genome.ucsc.edu/goldenPath/help/api.html).
#' - [RefSeq genomes](https://ftp.ncbi.nlm.nih.gov/genomes/refseq/).
#'
#' @examples
#' currentEnsemblBuild("Homo sapiens")
#' currentGencodeBuild("Homo sapiens")
#' currentRefSeqGenomeBuild(
#'     organism = "Homo sapiens",
#'     subdir = "vertebrate_mammalian"
#' )
#' currentUCSCGenomeBuild("Homo sapiens")
NULL



#' Get JSON (from REST API)
#'
#' Consider exporting this in pipette in a future update.
#'
#' @note Updated 2020-01-05.
#' @noRd
.getJSON <- function(url) {
    reponse <- GET(url = url, content_type("application/json"))
    stop_for_status(reponse)
    json <- fromJSON(toJSON(content(reponse)))
    assert(is.list(json))
    json
}



## Updated 2020-01-05.
#' @rdname currentGenomeBuild
#' @export
currentEnsemblBuild <- function(organism) {
    assert(isString(organism))
    organism <- snakeCase(organism)
    json <- .getJSON(pasteURL(
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



## Updated 2020-01-05.
#' @rdname currentGenomeBuild
#' @export
currentGencodeBuild <- function(organism) {
    organism <- match.arg(
        arg = organism,
        choices = c("Homo sapiens", "Mus musculus")
    )
    currentEnsemblBuild(organism)
}



## Alternate approach using URL only:
## https://ftp.ncbi.nlm.nih.gov/genomes/refseq/<subdir>/<organism>/
##     latest_assembly_versions/

#' @rdname currentGenomeBuild
#' @export
currentRefSeqGenomeBuild <- function(
    organism,
    subdir = "vertebrate_mammalian"
) {
    assert(
        isString(organism),
        isString(subdir)
    )
    url <- pasteURL(
        "ftp.ncbi.nlm.nih.gov",
        "genomes",
        "refseq",
        snakeCase(subdir),
        gsub(pattern = " ", replacement = "_", x = organism),
        "assembly_summary.txt",
        protocol = "https"
    )
    df <- import(
        file = url,
        format = "tsv",
        colnames = FALSE,
        skip = 2L,
        quiet = TRUE
    )
    vec <- as.character(df[1L, , drop = TRUE])
    if (!isSubset(organism, vec)) {
        stop("Invalid organism.")
    }
    out <- vec[[1L]]
    assert(isString(out))
    out
}



## Updated 2020-01-05.
#' @rdname currentGenomeBuild
#' @export
currentUCSCGenomeBuild <- function(organism) {
    assert(isString(organism))
    url <- "https://api.genome.ucsc.edu/list/ucscGenomes"
    json <- .getJSON(url)
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
