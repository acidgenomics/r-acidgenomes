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
`updateObject,EnsemblGenes` <- # nolint
    function(object, ..., verbose = FALSE) {
        stop("FIXME")

        new(Class = "EnsemblGenes", object)
    }



#' @rdname updateObject
#' @export
setMethod(
    f = "updateObject",
    signature = signature(object = "EnsemblGenes"),
    definition = `updateObject,EnsemblGenes`
)
