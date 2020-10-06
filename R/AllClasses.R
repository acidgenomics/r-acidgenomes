#' Ensembl-to-Entrez gene identifier mappings
#'
#' @details
#' Contains a `DataFrame` with `ensembl` and `entrez` columns.
#'
#' @export
#' @note Updated 2020-10-05.
#'
#' @return `Ensembl2Entrez`.
setClass(
    Class = "Ensembl2Entrez",
    contains = "DataFrame"
)
setValidity(
    Class = "Ensembl2Entrez",
    method = function(object) {
        validate(
            hasRows(object),
            identical(c("ensembl", "entrez"), colnames(object)),
            is.integer(object[["entrez"]])
        )
    }
)



#' Entrez-to-Ensembl gene identifier mappings
#'
#' @inherit Ensembl2Entrez-class details
#'
#' @export
#' @note Updated 2020-10-05.
#'
#' @return `Entrez2Ensembl`.
setClass(
    Class = "Entrez2Ensembl",
    contains = "DataFrame"
)
setValidity(
    Class = "Entrez2Ensembl",
    method = function(object) {
        validate(
            hasRows(object),
            identical(c("entrez", "ensembl"), colnames(object)),
            is.integer(object[["entrez"]])
        )
    }
)



#' Gene-to-symbol mappings
#'
#' @details
#' Contains a `DataFrame` with `geneID` and `geneName` columns.
#'
#' @section Genome metadata:
#'
#' We recommend slotting `organism`, `genomeBuild`, and `ensemblRelease` into
#' [`metadata()`][S4Vectors::metadata].
#'
#' @export
#' @note Updated 2019-08-08.
#'
#' @return `Gene2Symbol`.
setClass(
    Class = "Gene2Symbol",
    contains = "DataFrame"
)
setValidity(
    Class = "Gene2Symbol",
    method = function(object) {
        validate(
            hasRows(object),
            identical(colnames(object), c("geneID", "geneName")),
            is.character(object[["geneID"]]),
            isAny(object[["geneName"]], c("character", "factor"))
        )
    }
)



#' HGNC complete set metadata
#'
#' @export
#' @note Updated 2020-10-05.
#'
#' @return `HGNC`.
setClass(
    Class = "HGNC",
    contains = "DataFrame"
)
setValidity(
    Class = "HGNC",
    method = function(object) {
        validate(
            hasRows(object),
            isSubset(c("hgncID", "ensemblGeneID"), colnames(object)),
            is.integer(object[["hgncID"]]),
            hasNoDuplicates(object[["hgncID"]])
        )
    }
)



#' HGNC-to-Ensembl gene identifier mappings
#'
#' @details
#' Contains a `DataFrame` with `hgnc` and `ensembl` columns.
#'
#' @export
#' @note Updated 2020-10-05.
#'
#' @return `HGNC2Ensembl`.
setClass(
    Class = "HGNC2Ensembl",
    contains = "DataFrame"
)
setValidity(
    Class = "HGNC2Ensembl",
    method = function(object) {
        validate(
            hasRows(object),
            identical(c("hgnc", "ensembl"), colnames(object)),
            is.integer(object[["hgnc"]]),
            hasNoDuplicates(object[["hgnc"]])
        )
    }
)



#' MGI-to-Ensembl gene identifier mappings
#'
#' @details
#' Contains a `DataFrame` with `mgi` and `ensembl` columns.
#'
#' @export
#' @note Updated 2020-10-05.
#'
#' @return `MGI2Ensembl`.
setClass(
    Class = "MGI2Ensembl",
    contains = "DataFrame"
)
setValidity(
    Class = "MGI2Ensembl",
    method = function(object) {
        validate(
            hasRows(object),
            identical(colnames(object), c("mgi", "ensembl")),
            is.integer(object[["mgi"]]),
            hasNoDuplicates(object[["mgi"]])
        )
    }
)



#' Protein-to-gene mappings
#'
#' @details
#' Contains a `DataFrame` with `proteinID`, `geneID`, and `geneName` columns.
#'
#' @section Genome metadata:
#'
#' We recommend slotting `organism`, `genomeBuild`, and `ensemblRelease` into
#' [`metadata()`][S4Vectors::metadata].
#'
#' @export
#' @note Updated 2020-10-05.
#'
#' @return `Protein2Gene`.
setClass(
    Class = "Protein2Gene",
    contains = "DataFrame"
)
setValidity(
    Class = "Protein2Gene",
    method = function(object) {
        validate(
            hasRows(object),
            identical(
                x = colnames(object),
                y = c("proteinID", "geneID", "geneName")
            ),
            all(bapply(X = object, FUN = is.character))
        )
    }
)



#' Transcript-to-gene identifier mappings
#'
#' @details
#' Contains a `DataFrame` with `transcriptID` and `geneID` columns.
#'
#' @section Genome metadata:
#'
#' We recommend slotting `organism`, `genomeBuild`, and `ensemblRelease` into
#' `metadata`.
#'
#' @export
#' @note Updated 2020-10-05.
#'
#' @return `Tx2Gene`.
setClass(
    Class = "Tx2Gene",
    contains = "DataFrame"
)
setValidity(
    Class = "Tx2Gene",
    method = function(object) {
        validate(
            hasRows(object),
            identical(colnames(object), c("transcriptID", "geneID")),
            all(vapply(
                X = object,
                FUN = is.character,
                FUN.VALUE = logical(1L)
            )),
            hasNoDuplicates(object[["transcriptID"]])
        )
    }
)
