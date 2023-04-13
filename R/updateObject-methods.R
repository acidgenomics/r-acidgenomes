## FIXME Rename `entrezId` to `ncbiGeneId`.
## FIXME Add updateObject support for EnsemblGenes.
## FIXME Also add updateObject support for EnsemblTranscripts here.
## FIXME Add update support for GencodeGenes and GencodeTranscripts.



#' Update object
#'
#' @name updateObject
#' @author Michael Steinbaugh
#' @note Updated 2023-03-07.
NULL



## Updated 2023-04-13.
.updateGRanges <-
    function(object, ..., verbose = FALSE) {
        if (isSubset("entrezId", colnames(mcols(object)))) {
            if (isTRUE(verbose)) {
                alert(sprintf(
                    "Renaming {.var %s} to {.var %s}.",
                    "entrezId", "ncbiGeneId"
                ))
            }
            mcols(object)[["ncbiGeneId"]] <- mcols(object)[["entrezId"]]
            mcols(object)[["entrezId"]] <- NULL
        }
        validObject(object)
        object
    }



## Updated 2023-04-13.
`updateObject,EnsemblGenes` <- # nolint
    .updateGRanges

## Updated 2023-04-13.
`updateObject,EnsemblTranscripts` <- # nolint
    .updateGRanges

## Updated 2023-04-13.
`updateObject,GencodeGenes` <- # nolint
    .updateGRanges

## Updated 2023-04-13.
`updateObject,GencodeTranscripts` <- # nolint
    .updateGRanges



#' @rdname updateObject
#' @export
setMethod(
    f = "updateObject",
    signature = signature(object = "EnsemblGenes"),
    definition = `updateObject,EnsemblGenes`
)

#' @rdname updateObject
#' @export
setMethod(
    f = "updateObject",
    signature = signature(object = "EnsemblTranscripts"),
    definition = `updateObject,EnsemblTranscripts`
)

#' @rdname updateObject
#' @export
setMethod(
    f = "updateObject",
    signature = signature(object = "GencodeGenes"),
    definition = `updateObject,GencodeGenes`
)

#' @rdname updateObject
#' @export
setMethod(
    f = "updateObject",
    signature = signature(object = "GencodeTranscripts"),
    definition = `updateObject,GencodeTranscripts`
)
