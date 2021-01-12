## FIXME RETHINK APPROACH WITH THESE, IMPROVING SUPPORT FOR REFSEQ.



#' Prototype metadata
#'
#' @note Updated 2020-10-06.
#' @noRd
.prototypeMetadata <- list(
    version = .version,
    date = Sys.Date()
)



## FIXME RETHINK `ensemblRelease` usage here

#' Slot genome metadata
#'
#' @note Updated 2020-10-06.
#' @noRd
.slotGenomeMetadata <- function(object) {
    metadata <- metadata(object)
    proto <- c(
        .prototypeMetadata,
        list(
            organism = character(),
            genomeBuild = character(),
            ensemblRelease = integer()
        )
    )
    proto <- proto[setdiff(names(proto), names(metadata))]
    c(proto, metadata)
}



#' Slot organism into metadata, if necessary
#' @note Updated 2020-10-06.
#' @noRd
.slotOrganism <- function(object) {
    if (is.null(metadata(object)[["organism"]])) {
        metadata(object)[["organism"]] <- tryCatch(
            expr = detectOrganism(object),
            error = function(e) character()
        )
    }
    object
}
