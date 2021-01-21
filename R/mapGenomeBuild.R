#' Interconvert between EBI (Ensembl), NCBI, and UCSC genome build names
#'
#' @name mapGenomeBuild
#' @note Updated 2021-01-20.
#'
#' @inheritParams AcidRoxygen::params
#'
#' @return `character`. Ensembl genome build as the value, UCSC build as the
#'   name. Stops on match failure.
#'
#' @seealso
#' - [UCSC hgGateway](https://genome.ucsc.edu/cgi-bin/hgGateway).
#' - [bcbio genome recipes](https://goo.gl/164J2P).
#'
#' @examples
#' from <- c("hg19", "hg38")
#' to <- mapUCSCBuildToEnsembl(from)
#' print(to)
NULL



## Updated 2021-01-20.
.mapGenomeBuild <- function(object, from, to) {
    assert(
        isCharacter(object),
        isString(from),
        isString(to)
    )
    map <- import(
        file = system.file(
            "extdata", "map-genome-build.rds",
            package = packageName()
        ),
        quiet = TRUE
    )
    match <- match(x = object, table = map[[tolower(from)]])
    if (any(is.na(match))) {
        stop(sprintf(
            "Failed to match %s genome build to %s: %s.",
            from, to,
            toString(object[which(is.na(match))], width = 100L)
        ))
    }
    out <- map[[tolower(to)]][match]
    names(out) <- map[[tolower(from)]][match]
    out
}



## Updated 2021-01-20.
mapNCBIBuildToUCSC <- function(object) {
    .mapGenomeBuild(object, from = "NCBI", to = "UCSC")
}



## Updated 2021-01-20.
mapUCSCBuildToNCBI <- function(object) {
    .mapGenomeBuild(object, from = "UCSC", to = "NCBI")
}



#' @rdname mapGenomeBuild
#' @export
mapUCSCBuildToEnsembl <- function(...) {
    mapUCSCBuildToNCBI(...)
}



#' @rdname mapGenomeBuild
#' @export
mapEnsemblBuildToUCSC <- function(...) {
    mapNCBIBuildToUCSC(...)
}
