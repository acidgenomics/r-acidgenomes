## FIXME REWORK THIS FILE.



## Unordered list.
.li <- "  -"



.prototypeMetadata <- list(
    version = packageVersion("basejump"),
    date = Sys.Date()
)



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



## Slot organism into `metadata()`, if necessary.
.slotOrganism <- function(object) {
    if (is.null(metadata(object)[["organism"]])) {
        metadata(object)[["organism"]] <- tryCatch(
            expr = detectOrganism(object),
            error = function(e) character()
        )
    }
    object
}
