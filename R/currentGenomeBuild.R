#' Current genome build
#'
#' @name currentGenomeBuild
#' @note Updated 2020-01-05.
#'
#' @inheritParams AcidRoxygen::params
#'
#' @return `character(1)`.
#'   Genome build.
#'
#' @seealso
#' - [Ensembl REST API](https://rest.ensembl.org).
#' - [UCSC Genome Browser REST API](https://genome.ucsc.edu/goldenPath/help/api.html).
#'
#' @examples
#' currentGencodeBuild("Homo sapiens")
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
    df <- purrr::map_df(.x = l, .f = unlist)
    df <- df[df[["active"]] == 1L, ]
    if (!isSubset(organism, unique(df[["scientificName"]]))) {
        stop("Invalid organism.")
    }
    df <- df[df[["scientificName"]] == organism, ]
    ## The latest genome build has the lowest "orderKey" value.
    df <- df[order(df[["orderKey"]]), ]
    out <- df[["build"]][[1L]]
    assert(isString(out))
    out
}



## FIXME WORMBASE REST API.
## Updated 2020-01-05.
#' @rdname currentGenomeBuild
#' @export
currentWormBaseGenomeBuild <- function(organism) {
    print("FIXME")
}



## FIXME FLYBASE REST API.
## FIXME WORMBASE REST API.
## Updated 2020-01-05.
#' @rdname currentGenomeBuild
#' @export
currentFlyBaseGenomeBuild <- function(organism) {
    print("FIXME")
}
