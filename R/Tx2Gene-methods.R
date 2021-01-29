#' @name Tx2Gene
#' @inherit Tx2Gene-class title description return
#'
#' @note No attempt is made to arrange the rows by transcript identifier.
#' @note Updated 2021-01-29.
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



## Updated 2020-01-29.
`Tx2Gene,data.frame` <-  # nolint
    function(object) {
        assert(identical(ncol(object), 2L))
        colnames(object) <- c("txId", "geneId")
        object <- as(object, "DataFrame")
        Tx2Gene(object)
    }



#' @rdname Tx2Gene
#' @export
setMethod(
    f = "Tx2Gene",
    signature = signature("data.frame"),
    definition = `Tx2Gene,data.frame`
)



## Updated 2021-01-29.
`Tx2Gene,DataFrame` <-  # nolint
    function(object) {
        meta <- metadata(object)
        colnames(object) <- camelCase(colnames(object), strict = TRUE)
        colnames(object) <- gsub(
            pattern = "^transcript",
            replacement = "tx",
            x = colnames(object)
        )
        cols <- c("txId", "geneId")
        assert(
            isSubset(cols, colnames(object)),
            hasRows(object)
        )
        df <- object[, cols, drop = FALSE]
        ## This ensures any `CharacterList` columns from `GRangesList` get
        ## coerced to `character` correctly.
        df[[1L]] <- as.character(df[[1L]])
        df[[2L]] <- as.character(df[[2L]])
        ## This step is needed for handling raw GFF annotations.
        df <- unique(df)
        ## Ensure identifiers return sorted, which is more intuitive.
        idx <- order(df)
        df <- df[idx, , drop = FALSE]
        metadata(df) <- meta
        new(Class = "Tx2Gene", df)
    }



#' @rdname Tx2Gene
#' @export
setMethod(
    f = "Tx2Gene",
    signature = signature("DataFrame"),
    definition = `Tx2Gene,DataFrame`
)



## Updated 2021-01-29.
`Tx2Gene,GRanges` <-  # nolint
    function(object) {
        df <- as(object, "DataFrame")
        metadata(df) <- metadata(object)
        Tx2Gene(df)
    }



#' @rdname Tx2Gene
#' @export
setMethod(
    f = "Tx2Gene",
    signature = signature("GRanges"),
    definition = `Tx2Gene,GRanges`
)



## FIXME THIS CURRENTLY DOESNT WORK FOR REFSEQ, BECAUSE TX AND GENE IDENTIFERS
## ARE THE SAME HERE...MESSED UP.

## Updated 2021-01-29.
`Tx2Gene,CompressedGRangesList` <-  # nolint
    function(object) {
        x <- as(object, "CompressedGRangesList")
        x <- unname(x)
        gr <- unlist(x)
        assert(is(gr, "GRanges"))
        Tx2Gene(object = gr)
    }



#' @rdname Tx2Gene
#' @export
setMethod(
    f = "Tx2Gene",
    signature = signature("CompressedGRangesList"),
    definition = `Tx2Gene,CompressedGRangesList`
)



## Updated 2021-01-29.
`Tx2Gene,SummarizedExperiment` <-  # nolint
    function(object) {
        object <- rowData(object, use.names = TRUE)
        Tx2Gene(object)
    }



#' @rdname Tx2Gene
#' @export
setMethod(
    f = "Tx2Gene",
    signature = signature("SummarizedExperiment"),
    definition = `Tx2Gene,SummarizedExperiment`
)
