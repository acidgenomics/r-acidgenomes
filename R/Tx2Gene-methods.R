#' @inherit Tx2Gene-class title description return
#' @name Tx2Gene
#'
#' @note No attempt is made to arrange the rows by transcript identifier.
#' @note Updated 2020-01-07.
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



## Updated 2020-01-07.
`Tx2Gene,matrix` <-  # nolint
    function(object) {
        assert(
            is.character(object),
            identical(ncol(object), 2L)
        )
        object <- as.data.frame(object, stringsAsFactors = FALSE)
        Tx2Gene(object)
    }



#' @rdname Tx2Gene
#' @export
setMethod(
    f = "Tx2Gene",
    signature = signature("matrix"),
    definition = `Tx2Gene,matrix`
)



## Updated 2020-01-07.
`Tx2Gene,data.frame` <-  # nolint
    function(object) {
        assert(identical(ncol(object), 2L))
        colnames(object) <- c("transcriptID", "geneID")
        object <- as(object, "DataFrame")
        Tx2Gene(object = object, metadata = FALSE)
    }



#' @rdname Tx2Gene
#' @export
setMethod(
    f = "Tx2Gene",
    signature = signature("data.frame"),
    definition = `Tx2Gene,data.frame`
)



## Updated 2020-01-07.
`Tx2Gene,DataFrame` <-  # nolint
    function(object, metadata = TRUE) {
        assert(isFlag(metadata))
        cols <- c("transcriptID", "geneID")
        assert(
            isSubset(cols, colnames(object)),
            hasRows(object)
        )
        df <- object[, cols, drop = FALSE]
        df <- decode(df)
        ## Ensure identifiers return sorted, which is more intuitive.
        idx <- order(df)
        df <- df[idx, , drop = FALSE]
        if (isTRUE(metadata)) {
            metadata(df) <- .slotGenomeMetadata(object)
        }
        new(Class = "Tx2Gene", df)
    }



#' @rdname Tx2Gene
#' @export
setMethod(
    f = "Tx2Gene",
    signature = signature("DataFrame"),
    definition = `Tx2Gene,DataFrame`
)



## Updated 2021-01-07.
`Tx2Gene,GRanges` <-  # nolint
    function(object) {
        df <- as(object, "DataFrame")
        ## This step is needed for handling raw GFF annotations.
        df <- unique(df)
        metadata(df) <- metadata(object)
        Tx2Gene(object = df, metadata = TRUE)
    }



#' @rdname Tx2Gene
#' @export
setMethod(
    f = "Tx2Gene",
    signature = signature("GRanges"),
    definition = `Tx2Gene,GRanges`
)



## Updated 2021-01-07.
`Tx2Gene,SummarizedExperiment` <-  # nolint
    function(object) {
        object <- rowData(object, use.names = TRUE)
        Tx2Gene(object = object, metadata = TRUE)
    }



#' @rdname Tx2Gene
#' @export
setMethod(
    f = "Tx2Gene",
    signature = signature("SummarizedExperiment"),
    definition = `Tx2Gene,SummarizedExperiment`
)
