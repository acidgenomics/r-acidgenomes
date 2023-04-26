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
}



`show,EnsemblGenes` <- # nolint
    .showGRanges

`show,EnsemblTranscripts` <- # nolint
    .showGRanges

`show,FlyBaseGenes` <- # nolint
    .showGRanges

`show,FlyBaseTranscripts` <- # nolint
    .showGRanges

`show,GencodeGenes` <- # nolint
    .showGRanges

`show,GencodeTranscripts` <- # nolint
    .showGRanges

`show,RefSeqGenes` <- # nolint
    .showGRanges

`show,RefSeqTranscripts` <- # nolint
    .showGRanges

`show,UCSCGenes` <- # nolint
    .showGRanges

`show,UCSCTranscripts` <- # nolint
    .showGRanges

`show,WormBaseGenes` <- # nolint
    .showGRanges

`show,WormBaseTranscripts` <- # nolint
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
    signature = signature(object = "FlyBaseGenes"),
    definition = `show,FlyBaseGenes`
)

#' @rdname show
#' @export
setMethod(
    f = "show",
    signature = signature(object = "FlyBaseTranscripts"),
    definition = `show,FlyBaseTranscripts`
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
    signature = signature(object = "RefSeqGenes"),
    definition = `show,RefSeqGenes`
)

#' @rdname show
#' @export
setMethod(
    f = "show",
    signature = signature(object = "RefSeqTranscripts"),
    definition = `show,RefSeqTranscripts`
)

#' @rdname show
#' @export
setMethod(
    f = "show",
    signature = signature(object = "UCSCGenes"),
    definition = `show,UCSCGenes`
)

#' @rdname show
#' @export
setMethod(
    f = "show",
    signature = signature(object = "UCSCTranscripts"),
    definition = `show,UCSCTranscripts`
)

#' @rdname show
#' @export
setMethod(
    f = "show",
    signature = signature(object = "WormBaseGenes"),
    definition = `show,WormBaseGenes`
)

#' @rdname show
#' @export
setMethod(
    f = "show",
    signature = signature(object = "WormBaseTranscripts"),
    definition = `show,WormBaseTranscripts`
)
