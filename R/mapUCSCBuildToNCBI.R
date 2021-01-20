#' Map UCSC genome build to NCBI
#'
#' @export
#' @note Updated 2021-01-20.
#'
#' @inheritParams AcidRoxygen::params
#'
#' @return `character`. Ensembl genome build as the value, UCSC build as the
#'   name. Stops on match failure.
#'
#' @seealso
#' - [UCSC hgGateway](https://genome.ucsc.edu/cgi-bin/hgGateway)
#' - [bcbio genome recipes](https://goo.gl/164J2P)
#'
#' @examples
#' from <- c("hg19", "hg38")
#' to <- mapUCSCBuildToEnsembl(from)
#' print(to)
mapUCSCBuildToNCBI <- function(object) {
    assert(isCharacter(object))
    map <- import(
        file = system.file(
            "extdata", "ncbi-to-ucsc.rds",
            package = packageName()
        ),
        quiet = TRUE
    )
    match <- match(x = object, table = map[["ucsc"]])
    ## Stop on any match failure.
    if (any(is.na(match))) {
        stop(sprintf(
            "Failed to match UCSC to Ensembl: %s.",
            toString(object[which(is.na(match))], width = 100L)
        ))
    }
    out <- map[["ensembl"]][match]
    names(out) <- map[["ucsc"]][match]
    out
}



## NOTE This is safe to deprecate and remove once we update basejump.

#' @rdname mapUCSCBuildToNCBI
#' @export
mapUCSCBuildToEnsembl <- function(...) {
    mapUCSCBuildToNCBI(...)
}
