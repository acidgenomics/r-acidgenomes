#' Show an object
#'
#' @name show
#' @author Michael Steinbaugh
#' @note Updated 2023-04-26.
#'
#' @inheritParams AcidRoxygen::params
#'
#' @return Console output.
NULL



## Updated 2023-04-26.
.showGRanges <- function(object) {
    assert(validObject(object))
    showHeader(object)
    meta <- metadata(object)
    showSlotInfo(list(
        "organism" = meta[["organism"]],
        "genomeBuild" = meta[["genomeBuild"]],
        "release" = meta[["release"]],
        "ignoreVersion" = meta[["ignoreVersion"]],
        "names" = names(object)
    ))
}



`show,EnsemblGenes` <- # nolint
    .showGRanges

`show,EnsemblTranscripts` <- # nolint
    .showGRanges

`show,FlybaseGenes` <- # nolint
    .showGRanges

`show,FlybaseTranscripts` <- # nolint
    .showGRanges

`show,GencodeGenes` <- # nolint
    .showGRanges

`show,GencodeTranscripts` <- # nolint
    .showGRanges

`show,RefseqGenes` <- # nolint
    .showGRanges

`show,RefseqTranscripts` <- # nolint
    .showGRanges

`show,UcscGenes` <- # nolint
    .showGRanges

`show,UcscTranscripts` <- # nolint
    .showGRanges

`show,WormbaseGenes` <- # nolint
    .showGRanges

`show,WormbaseTranscripts` <- # nolint
    .showGRanges



#' @rdname show
#' @export
setMethod(
    f = "show",
    signature = signature(object = "EnsemblGenes"),
    definition = `show,EnsemblGenes`
)

#' @rdname show
#' @export
setMethod(
    f = "show",
    signature = signature(object = "EnsemblTranscripts"),
    definition = `show,EnsemblTranscripts`
)

#' @rdname show
#' @export
setMethod(
    f = "show",
    signature = signature(object = "FlybaseGenes"),
    definition = `show,FlybaseGenes`
)

#' @rdname show
#' @export
setMethod(
    f = "show",
    signature = signature(object = "FlybaseTranscripts"),
    definition = `show,FlybaseTranscripts`
)

#' @rdname show
#' @export
setMethod(
    f = "show",
    signature = signature(object = "GencodeGenes"),
    definition = `show,GencodeGenes`
)

#' @rdname show
#' @export
setMethod(
    f = "show",
    signature = signature(object = "GencodeTranscripts"),
    definition = `show,GencodeTranscripts`
)

#' @rdname show
#' @export
setMethod(
    f = "show",
    signature = signature(object = "RefseqGenes"),
    definition = `show,RefseqGenes`
)

#' @rdname show
#' @export
setMethod(
    f = "show",
    signature = signature(object = "RefseqTranscripts"),
    definition = `show,RefseqTranscripts`
)

#' @rdname show
#' @export
setMethod(
    f = "show",
    signature = signature(object = "UcscGenes"),
    definition = `show,UcscGenes`
)

#' @rdname show
#' @export
setMethod(
    f = "show",
    signature = signature(object = "UcscTranscripts"),
    definition = `show,UcscTranscripts`
)

#' @rdname show
#' @export
setMethod(
    f = "show",
    signature = signature(object = "WormbaseGenes"),
    definition = `show,WormbaseGenes`
)

#' @rdname show
#' @export
setMethod(
    f = "show",
    signature = signature(object = "WormbaseTranscripts"),
    definition = `show,WormbaseTranscripts`
)
