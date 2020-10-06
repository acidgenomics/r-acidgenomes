#' Map UCSC genome build to Ensembl
#'
#' @export
#' @note Updated 2020-10-06.
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
mapUCSCBuildToEnsembl <- function(object) {
    assert(isCharacter(object))
    map <- import(
        file = system.file(
            "extdata", "mapUCSCBuildToEnsembl.rds",
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
