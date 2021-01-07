## FIXME CONSIDER SORTING AND MAKING THESE UNIQUE AS AN AUTOMATIC STEP HERE.
## FIXME OBJECT CLASS NEEDS TO CHECK THAT THERE ARE NO DUPLICATES.
## FIXME CONSIDER UPDATING AND SORTING LEGACY OBJECTS AUTOMATICALLY?



#' @inherit Tx2Gene-class title description return
#' @name Tx2Gene
#'
#' @note No attempt is made to arrange the rows by transcript identifier.
#' @note Updated 2019-10-24.
#'
#' @inheritParams AcidRoxygen::params
#' @param metadata `logical(1)`.
#'   Include genome metadata.
#'
#' @seealso [makeTx2Gene()].
#'
#' @examples
#' ## SummarizedExperiment ====
#' data(SummarizedExperiment_transcripts, package = "AcidTest")
#' txse <- SummarizedExperiment_transcripts
#'
#' ## SummarizedExperiment ====
#' x <- Tx2Gene(txse)
#' print(x)
NULL



## Updated 2019-10-24.
`Tx2Gene,matrix` <-  # nolint
    function(object) {
        assert(
            is.character(object),
            identical(ncol(object), 2L)
        )
        Tx2Gene(object = as.data.frame(object, stringsAsFactors = FALSE))
    }



#' @rdname Tx2Gene
#' @export
setMethod(
    f = "Tx2Gene",
    signature = signature("matrix"),
    definition = `Tx2Gene,matrix`
)



## Updated 2019-10-24.
`Tx2Gene,data.frame` <-  # nolint
    function(object) {
        assert(identical(ncol(object), 2L))
        colnames(object) <- c("transcriptID", "geneID")
        Tx2Gene(object = as(object, "DataFrame"), metadata = FALSE)
    }



#' @rdname Tx2Gene
#' @export
setMethod(
    f = "Tx2Gene",
    signature = signature("data.frame"),
    definition = `Tx2Gene,data.frame`
)



## Updated 2019-10-24.
`Tx2Gene,DataFrame` <-  # nolint
    function(object, metadata = TRUE) {
        assert(isFlag(metadata))
        cols <- c("transcriptID", "geneID")
        assert(
            isSubset(cols, colnames(object)),
            hasRows(object)
        )
        data <- decode(object[, cols, drop = FALSE])
        if (isTRUE(metadata)) {
            metadata(data) <- .slotGenomeMetadata(object)
        }
        new(Class = "Tx2Gene", data)
    }



#' @rdname Tx2Gene
#' @export
setMethod(
    f = "Tx2Gene",
    signature = signature("DataFrame"),
    definition = `Tx2Gene,DataFrame`
)



## Updated 2019-10-24.
`Tx2Gene,GRanges` <-  # nolint
    function(object) {
        data <- as(object, "DataFrame")
        ## This step is needed for handling raw GFF annotations.
        data <- unique(data)
        metadata(data) <- metadata(object)
        Tx2Gene(object = data, metadata = TRUE)
    }



#' @rdname Tx2Gene
#' @export
setMethod(
    f = "Tx2Gene",
    signature = signature("GRanges"),
    definition = `Tx2Gene,GRanges`
)



## Updated 2019-10-24.
`Tx2Gene,SummarizedExperiment` <-  # nolint
    function(object) {
        Tx2Gene(object = rowData(object, use.names = TRUE), metadata = TRUE)
    }



#' @rdname Tx2Gene
#' @export
setMethod(
    f = "Tx2Gene",
    signature = signature("SummarizedExperiment"),
    definition = `Tx2Gene,SummarizedExperiment`
)
