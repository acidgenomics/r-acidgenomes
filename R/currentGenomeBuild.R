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
#' - [UCSC REST API](https://genome.ucsc.edu/goldenPath/help/api.html).
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



## Updated 2020-01-05.
#' @rdname currentGenomeBuild
#' @export
currentEnsemblBuild <- function(organism) {
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



## FIXME SUBDIR INPUT HERE IS ANNOYING. TAKE THIS STEP OUT.

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
        isOrganism(organism),
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



#' Get the RefSeq organism directory structure on the FTP server
#'
#' @note Updated 2021-01-08.
#' @noRd
.getRefSeqGenomeURL <- function(
    organism,
    subdir = NULL,
    quiet = FALSE
) {
    assert(
        isOrganism(organism),
        isString(subdir, nullOK = TRUE),
        isFlag(quiet)
    )
    if (isFALSE(quiet)) {
        alert(sprintf(
            "Locating {.emph %s} genome on RefSeq FTP server.",
            organism
        ))
    }
    baseURL <- "ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq"
    ## Match the organism to the top-level subdirectory, if necessary.
    if (is.null(subdir)) {
        subdirs <- getURLDirList(url = baseURL)
        keep <- grepl(pattern = "^[a-z_]+$", x = subdirs)
        subdirs <- sort(subdirs[keep])
        list <- bplapply(
            X = subdirs,
            baseURL = baseURL,
            FUN = function(subdir, baseURL) {
                url <- pasteURL(baseURL, subdir)
                x <- getURLDirList(url = url)
                keep <- grepl(pattern = "^[A-Z][a-z]+_[a-z]+$", x = x)
                x <- sort(x[keep])
                x
            }
        )
        names(list) <- subdirs
        match <- vapply(
            X = list,
            organism = gsub(pattern = " ", replacement = "_", x = organism),
            FUN = function(strings, organism) {
                isSubset(x = organism, y = strings)
            },
            FUN.VALUE = logical(1L),
            USE.NAMES = TRUE
        )
        subdir <- names(match)[match]
        assert(isString(subdir))
    }
    url <- pasteURL(
        baseURL, subdir, gsub(pattern = " ", replacement = "_", x = organism)
    )
    assert(url.exists(url))
    if (isFALSE(quiet)) {
        dl(c("URL" = url))
    }
    url
}
