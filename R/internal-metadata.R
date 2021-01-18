#' Prototype metadata
#'
#' @note Updated 2020-10-06.
#' @noRd
.prototypeMetadata <- list(
    "date" = Sys.Date(),
    "version" = .version
)



#' Slot genome metadata
#'
#' @note Updated 2021-01-18.
#' @noRd
.slotGenomeMetadata <- function(object) {
    metadata <- metadata(object)
    proto <- append(
        x = .prototypeMetadata,
        values = list(
            "genomeBuild" = character(),
            "organism" = character(),
            "release" = integer()
        )
    )
    proto <- proto[setdiff(names(proto), names(metadata))]
    out <- append(x = proto, values = metadata)
    out <- out[sort(unique(names(out)))]
    out
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
