#' Organism
#'
#' @name organism
#' @note Updated 2021-02-01.
#'
#' @inheritParams AcidRoxygen::params
#'
#' @return `character(1)`.
#' Latin organism name (e.g. *Homo sapiens*).
#'
#' @seealso [detectOrganism()][.
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



#' @rdname organism
#' @export
setMethod(
    f = "organism",
    signature = signature("matrix"),
    definition = `organism,matrix`
)



## Updated 2020-01-30.
`organism,Matrix` <-  # nolint
    `organism,matrix`



#' @rdname organism
#' @export
setMethod(
    f = "organism",
    signature = signature("Matrix"),
    definition = `organism,Matrix`
)



## Updated 2019-07-22.
`organism,data.frame` <-  # nolint
    `organism,matrix`



#' @rdname organism
#' @export
setMethod(
    f = "organism",
    signature = signature("data.frame"),
    definition = `organism,data.frame`
)



## Note that DataFrame and GRanges inherit from this class.
## Updated 2019-07-22.
`organism,Annotated` <-  # nolint
    function(object) {
        metadata(object)[["organism"]]
    }



#' @rdname organism
#' @export
setMethod(
    f = "organism",
    signature = signature("Annotated"),
    definition = `organism,Annotated`
)



## Updated 2019-07-22.
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



#' @rdname organism
#' @export
setMethod(
    f = "organism",
    signature = signature("DataFrame"),
    definition = `organism,DataFrame`
)



## Updated 2021-02-02.
`organism,GRanges` <-  # nolint
    function(object) {
        ## Attempt to use metadata stash, if defined.
        organism <- `organism,Annotated`(object)
        if (isString(organism)) {
            return(organism)
        }
        assert(hasNames(object))
        detectOrganism(names(object))
    }



#' @rdname organism
#' @export
setMethod(
    f = "organism",
    signature = signature("GRanges"),
    definition = `organism,GRanges`
)



## Updated 2019-08-06.
`organism<-,Annotated,character` <-  # nolint
    function(object, value) {
        metadata(object)[["organism"]] <- value
        object
    }



#' @rdname organism
#' @export
setReplaceMethod(
    f = "organism",
    signature = signature("Annotated"),
    definition = `organism<-,Annotated,character`
)
