#' Organism
#'
#' @name organism
#' @note Updated 2021-08-03.
#'
#' @inheritParams AcidRoxygen::params
#'
#' @return `character(1)`.
#' Latin organism name (e.g. *Homo sapiens*).
#'
#' @seealso `detectOrganism()`.
#'
#' @examples
#' ## DataFrame ====
#' df <- DataFrame(
#'     "txId" = c(
#'         "ENST00000635602.1",
#'         "ENST00000635506.1"
#'     ),
#'     "geneId" = c(
#'         "ENSG00000283061.1",
#'         "ENSG00000283061.1"
#'     ),
#'     row.names = c(
#'         "ENST00000635602.1",
#'         "ENST00000635506.1"
#'     )
#' )
#' x <- organism(df)
#' print(x)
NULL



## Assuming gene identifiers are defined in the rownames.
## Updated 2019-07-22.
`organism,matrix` <-  # nolint
    function(object) {
        assert(hasRownames(object))
        detectOrganism(rownames(object))
    }



## Updated 2020-01-30.
`organism,Matrix` <-  # nolint
    `organism,matrix`



## Updated 2019-07-22.
`organism,data.frame` <-  # nolint
    `organism,matrix`



## Note that DataFrame and GenomicRanges inherit from this class.
## Updated 2019-07-22.
`organism,Annotated` <-  # nolint
    function(object) {
        metadata(object)[["organism"]]
    }



## Updated 2021-10-13.
`organism,DataFrame` <-  # nolint
    function(object) {
        ## Attempt to use metadata stash, if defined.
        organism <- `organism,Annotated`(object)
        if (isOrganism(organism)) {
            return(organism)
        }
        ## Otherwise, fall back to matrix method.
        `organism,matrix`(object)
    }



## Updated 2021-02-02.
`organism,GenomicRanges` <-  # nolint
    function(object) {
        ## Attempt to use metadata stash, if defined.
        organism <- `organism,Annotated`(object)
        if (isString(organism)) {
            return(organism)
        }
        assert(hasNames(object))
        detectOrganism(names(object))
    }



## Updated 2019-08-06.
`organism<-,Annotated,character` <-  # nolint
    function(object, value) {
        metadata(object)[["organism"]] <- value
        object
    }



#' @rdname organism
#' @export
setMethod(
    f = "organism",
    signature = signature(object = "Annotated"),
    definition = `organism,Annotated`
)

#' @rdname organism
#' @export
setMethod(
    f = "organism",
    signature = signature(object = "DataFrame"),
    definition = `organism,DataFrame`
)

#' @rdname organism
#' @export
setMethod(
    f = "organism",
    signature = signature(object = "GenomicRanges"),
    definition = `organism,GenomicRanges`
)

#' @rdname organism
#' @export
setMethod(
    f = "organism",
    signature = signature(object = "Matrix"),
    definition = `organism,Matrix`
)

#' @rdname organism
#' @export
setMethod(
    f = "organism",
    signature = signature(object = "data.frame"),
    definition = `organism,data.frame`
)

#' @rdname organism
#' @export
setMethod(
    f = "organism",
    signature = signature(object = "matrix"),
    definition = `organism,matrix`
)

## NOTE Can't define "value" in signature here.
#' @rdname organism
#' @export
setReplaceMethod(
    f = "organism",
    signature = signature(object = "Annotated"),
    definition = `organism<-,Annotated,character`
)
