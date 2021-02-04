#' @rdname Ensembl2Entrez
#' @export
setGeneric(
    name = "Ensembl2Entrez",
    def = function(object, ...) {
        standardGeneric("Ensembl2Entrez")
    }
)



#' @rdname Ensembl2Entrez
#' @export
setGeneric(
    name = "Entrez2Ensembl",
    def = function(object, ...) {
        standardGeneric("Entrez2Ensembl")
    }
)



#' @rdname Gene2Symbol
#' @export
setGeneric(
    name = "Gene2Symbol",
    def = function(object, ...) {
        standardGeneric("Gene2Symbol")
    }
)



#' @rdname Tx2Gene
#' @export
setGeneric(
    name = "Tx2Gene",
    def = function(object, ...) {
        standardGeneric("Tx2Gene")
    }
)



#' @rdname organism
#' @name organism
#' @importFrom AcidGenerics organism
#' @usage organism(object)
#' @export
NULL

#' @rdname organism
#' @name organism<-
#' @importFrom AcidGenerics organism<-
#' @usage organism(object) <- value
#' @export
NULL



#' @rdname stripVersions
#' @name stripGeneVersions
#' @importFrom AcidGenerics stripGeneVersions
#' @usage stripGeneVersions(object, ...)
#' @export
NULL

#' @rdname stripVersions
#' @name stripTranscriptVersions
#' @importFrom AcidGenerics stripTranscriptVersions
#' @usage stripTranscriptVersions(object, ...)
#' @export
NULL



#' @rdname summary
#' @name summary
#' @importFrom AcidGenerics summary
#' @usage summary(object, ...)
#' @export
NULL
