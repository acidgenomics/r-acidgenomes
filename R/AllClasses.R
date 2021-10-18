## Genome annotation classes ===================================================

#' Shared GenomicRanges validity checks
#'
#' @note Note that genome build and organism are not defined in minimal FlyBase
#'   GTF example.
#' @note Updated 2021-03-03.
#' @noRd
.grangesValidity <- function(object) {
    if (is(object, "GenomicRangesList")) {
        gr <- object[[1L]]
    } else {
        gr <- object
    }
    ok <- validateClasses(
        object = metadata(object),
        expected = list(
            "date" = "Date",
            "ignoreVersion" = "logical",
            "level" = "character",
            "packageVersion" = "package_version",
            "provider" = "character",
            "synonyms" = "logical"
        ),
        subset = TRUE
    )
    if (!isTRUE(ok)) {
        return(ok)
    }
    TRUE
}



#' Shared Ensembl validity checks
#'
#' @note Updated 2021-02-26.
#' @noRd
.ensemblValidity <- function(object) {
    ok <- .grangesValidity(object)
    if (!isTRUE(ok)) {
        return(ok)
    }
    ok <- validate(
        identical(metadata(object)[["provider"]], "Ensembl")
    )
    if (!isTRUE(ok)) {
        return(ok)
    }
    ok <- validateClasses(
        object = metadata(object),
        expected = list(
            ## > "genomeBuild" = "character",
            "organism" = "character"
        ),
        subset = TRUE
    )
    if (!isTRUE(ok)) {
        return(ok)
    }
    ## Can't always get the release version from GFF file, so just check
    ## for ensembldb return.
    if (isSubset("ensembldb", names(metadata(object)))) {
        ok <- validate(
            is.integer(metadata(object)[["release"]])
        )
        if (!isTRUE(ok)) {
            return(ok)
        }
    }
    TRUE
}



#' Shared FlyBase validity checks
#'
#' @note Updated 2021-01-30.
#' @noRd
.flybaseValidity <- function(object) {
    ok <- .grangesValidity(object)
    if (!isTRUE(ok)) {
        return(ok)
    }
    ok <- validate(
        identical(metadata(object)[["provider"]], "FlyBase")
    )
    if (!isTRUE(ok)) {
        return(ok)
    }
    TRUE
}



#' Shared GENCODE validity checks
#'
#' @note Updated 2021-02-26.
#' @noRd
.gencodeValidity <- function(object) {
    ok <- .grangesValidity(object)
    if (!isTRUE(ok)) {
        return(ok)
    }
    ok <- validate(
        identical(metadata(object)[["provider"]], "GENCODE")
    )
    if (!isTRUE(ok)) {
        return(ok)
    }
    ok <- validateClasses(
        object = metadata(object),
        expected = list(
            "genomeBuild" = "character",
            "organism" = "character"
        ),
        subset = TRUE
    )
    if (!isTRUE(ok)) {
        return(ok)
    }
    TRUE
}



#' Shared RefSeq validity checks
#'
#' @note Updated 2021-02-01.
#' @noRd
.refseqValidity <- function(object) {
    ok <- .grangesValidity(object)
    if (!isTRUE(ok)) {
        return(ok)
    }
    ok <- validate(
        identical(metadata(object)[["provider"]], "RefSeq")
    )
    if (!isTRUE(ok)) {
        return(ok)
    }
    ok <- validateClasses(
        object = metadata(object),
        expected = list(
            "genomeBuild" = "character",
            "organism" = "character"
        ),
        subset = TRUE
    )
    if (!isTRUE(ok)) {
        return(ok)
    }
    TRUE
}



#' Shared UCSC validity checks
#'
#' @note Updated 2021-02-01.
#' @noRd
.ucscValidity <- function(object) {
    ok <- .grangesValidity(object)
    if (!isTRUE(ok)) {
        return(ok)
    }
    ok <- validate(
        identical(metadata(object)[["provider"]], "UCSC")
    )
    if (!isTRUE(ok)) {
        return(ok)
    }
    ok <- validateClasses(
        object = metadata(object),
        expected = list(
            "genomeBuild" = "character",
            "organism" = "character"
        ),
        subset = TRUE
    )
    if (!isTRUE(ok)) {
        return(ok)
    }
    TRUE
}



#' Shared WormBase validity checks
#'
#' @note Updated 2021-02-01.
#' @noRd
.wormbaseValidity <- function(object) {
    ok <- .grangesValidity(object)
    if (!isTRUE(ok)) {
        return(ok)
    }
    ok <- validate(
        identical(metadata(object)[["provider"]], "WormBase")
    )
    if (!isTRUE(ok)) {
        return(ok)
    }
    TRUE
}



#' Ensembl gene annotations
#'
#' @details
#' Contains `GenomicRanges` with Ensembl gene-level annotations.
#'
#' @export
#' @note Updated 2021-10-14.
#'
#' @return `EnsemblGenes`.
setClass(
    Class = "EnsemblGenes",
    contains = "GRanges"
)
setValidity(
    Class = "EnsemblGenes",
    method = function(object) {
        ok <- .ensemblValidity(object)
        if (!isTRUE(ok)) {
            return(ok)
        }
        ok <- validate(identical(metadata(object)[["level"]], "genes"))
        if (!isTRUE(ok)) {
            return(ok)
        }
        TRUE
    }
)



#' Ensembl transcript annotations
#'
#' @details
#' Contains `GenomicRanges` with Ensembl transcript-level annotations.
#'
#' @export
#' @note Updated 2021-10-14.
#'
#' @return `EnsemblTranscripts`.
setClass(
    Class = "EnsemblTranscripts",
    contains = "GRanges"
)
setValidity(
    Class = "EnsemblTranscripts",
    method = function(object) {
        ok <- .ensemblValidity(object)
        if (!isTRUE(ok)) {
            return(ok)
        }
        ok <- validate(identical(metadata(object)[["level"]], "transcripts"))
        if (!isTRUE(ok)) {
            return(ok)
        }
        TRUE
    }
)



#' FlyBase gene annotations
#'
#' @details
#' Contains `GenomicRanges` with FlyBase gene-level annotations.
#'
#' @export
#' @note Updated 2021-10-14.
#'
#' @return `FlyBaseGenes`.
setClass(
    Class = "FlyBaseGenes",
    contains = "GRanges"
)
setValidity(
    Class = "FlyBaseGenes",
    method = function(object) {
        ok <- .flybaseValidity(object)
        if (!isTRUE(ok)) {
            return(ok)
        }
        TRUE
    }
)



#' FlyBase transcript annotations
#'
#' @details
#' Contains `GenomicRanges` with FlyBase transcript-level annotations.
#'
#' @export
#' @note Updated 2021-10-14.
#'
#' @return `FlyBaseTranscripts`.
setClass(
    Class = "FlyBaseTranscripts",
    contains = "GRanges"
)
setValidity(
    Class = "FlyBaseTranscripts",
    method = function(object) {
        ok <- .flybaseValidity(object)
        if (!isTRUE(ok)) {
            return(ok)
        }
        TRUE
    }
)



#' GENCODE gene annotations
#'
#' @details
#' Contains `GenomicRanges` with GENCODE gene-level annotations.
#'
#' @export
#' @note Updated 2021-10-14.
#'
#' @return `GencodeGenes`.
setClass(
    Class = "GencodeGenes",
    contains = "GRanges"
)
setValidity(
    Class = "GencodeGenes",
    method = function(object) {
        ok <- .gencodeValidity(object)
        if (!isTRUE(ok)) {
            return(ok)
        }
        TRUE
    }
)



#' GENCODE transcript annotations
#'
#' @details
#' Contains `GenomicRanges` with GENCODE transcript-level annotations.
#'
#' @export
#' @note Updated 2021-10-14.
#'
#' @return `GencodeTranscripts`.
setClass(
    Class = "GencodeTranscripts",
    contains = "GRanges"
)
setValidity(
    Class = "GencodeTranscripts",
    method = function(object) {
        ok <- .gencodeValidity(object)
        if (!isTRUE(ok)) {
            return(ok)
        }
        TRUE
    }
)



#' RefSeq gene annotations
#'
#' @details
#' Contains a `GenomicRangesList` with RefSeq gene-level annotations.
#'
#' @export
#' @note Updated 2021-10-14.
#'
#' @return `RefSeqGenes`.
setClass(
    Class = "RefSeqGenes",
    contains = "CompressedGRangesList"
)
setValidity(
    Class = "RefSeqGenes",
    method = function(object) {
        ok <- .refseqValidity(object)
        if (!isTRUE(ok)) {
            return(ok)
        }
        TRUE
    }
)



#' RefSeq transcript annotations
#'
#' @details
#' Contains a `GenomicRangesList` with RefSeq transcript-level annotations.
#'
#' @export
#' @note Updated 2021-10-14.
#'
#' @return `RefSeqTranscripts`.
setClass(
    Class = "RefSeqTranscripts",
    contains = "CompressedGRangesList"
)
setValidity(
    Class = "RefSeqTranscripts",
    method = function(object) {
        ok <- .refseqValidity(object)
        if (!isTRUE(ok)) {
            return(ok)
        }
        TRUE
    }
)



#' UCSC gene annotations
#'
#' @details
#' Contains `GenomicRanges` with UCSC gene-level annotations.
#'
#' @export
#' @note Updated 2021-10-14.
#'
#' @return `UCSCGenes`.
setClass(
    Class = "UCSCGenes",
    contains = "GRanges"
)
setValidity(
    Class = "UCSCGenes",
    method = function(object) {
        ok <- .ucscValidity(object)
        if (!isTRUE(ok)) {
            return(ok)
        }
        TRUE
    }
)



#' UCSC transcript annotations
#'
#' @details
#' Contains `GenomicRanges` with UCSC transcript-level annotations.
#'
#' @export
#' @note Updated 2021-10-14.
#'
#' @return `UCSCTranscripts`.
setClass(
    Class = "UCSCTranscripts",
    contains = "GRanges"
)
setValidity(
    Class = "UCSCTranscripts",
    method = function(object) {
        ok <- .ucscValidity(object)
        if (!isTRUE(ok)) {
            return(ok)
        }
        TRUE
    }
)



#' WormBase gene annotations
#'
#' @details
#' Contains `GenomicRanges` with WormBase gene-level annotations.
#'
#' @export
#' @note Updated 2021-10-14.
#'
#' @return `WormBaseGenes`.
setClass(
    Class = "WormBaseGenes",
    contains = "GRanges"
)
setValidity(
    Class = "WormBaseGenes",
    method = function(object) {
        ok <- .wormbaseValidity(object)
        if (!isTRUE(ok)) {
            return(ok)
        }
        TRUE
    }
)



#' WormBase transcript annotations
#'
#' @details
#' Contains `GenomicRanges` with WormBase transcript-level annotations.
#'
#' @export
#' @note Updated 2021-10-14.
#'
#' @return `WormBaseTranscripts`.
setClass(
    Class = "WormBaseTranscripts",
    contains = "GRanges"
)
setValidity(
    Class = "WormBaseTranscripts",
    method = function(object) {
        ok <- .wormbaseValidity(object)
        if (!isTRUE(ok)) {
            return(ok)
        }
        TRUE
    }
)



## Identifier classes ==========================================================

#' NCBI Entrez gene identifier information
#'
#' @export
#' @note Updated 2021-02-12.
#'
#' @return `EntrezGeneInfo`.
setClass(
    Class = "EntrezGeneInfo",
    contains = "DataFrame"
)
setValidity(
    Class = "EntrezGeneInfo",
    method = function(object) {
        TRUE
    }
)



#' HGNC complete set metadata
#'
#' @export
#' @note Updated 2021-02-07.
#'
#' @return `HGNC`.
setClass(
    Class = "HGNC",
    contains = "DataFrame"
)
setValidity(
    Class = "HGNC",
    method = function(object) {
        ok <- validate(
            hasColnames(object),
            hasRows(object)
        )
        if (!isTRUE(ok)) {
            return(ok)
        }
        cols <- c(
            "hgnc" = "hgncId",
            "ensembl" = "ensemblGeneId"
        )
        if (!isSubset(cols, colnames(object))) {
            colnames(object) <- camelCase(colnames(object), strict = TRUE)
        }
        ok <- validate(
            isSubset(cols, colnames(object)),
            is.integer(object[[cols[[1L]]]]),
            hasNoDuplicates(object[[cols[[1L]]]])
        )
        if (!isTRUE(ok)) {
            return(ok)
        }
        TRUE
    }
)



## Identifier mapping classes ==================================================

#' @inherit AcidGenerics::Ensembl2Entrez description return title
#' @note Updated 2021-08-18.
#' @export
#'
#' @details
#' Contains a `DataFrame` with `ensemblId` and `entrezId` columns.
#'
#' @return `Ensembl2Entrez`.
setClass(
    Class = "Ensembl2Entrez",
    contains = "DataFrame"
)
setValidity(
    Class = "Ensembl2Entrez",
    method = function(object) {
        ok <- validate(
            identical(ncol(object), 2L),
            hasColnames(object),
            hasRows(object),
            all(complete.cases(object))
        )
        if (!isTRUE(ok)) {
            return(ok)
        }
        cols <- c(
            "ensembl" = "ensemblId",
            "entrez" = "entrezId"
        )
        if (!identical(cols, colnames(object))) {
            colnames(object) <- unname(cols)
        }
        ok <- validate(
            is.integer(object[[cols[["entrez"]]]])
        )
        if (!isTRUE(ok)) {
            return(ok)
        }
        TRUE
    }
)



## FIXME Inherit the documentation from AcidGenerics.

#' Entrez-to-Ensembl gene identifier mappings
#'
#' @inherit Ensembl2Entrez-class details
#'
#' @export
#' @note Updated 2021-08-10.
#'
#' @return `Entrez2Ensembl`.
setClass(
    Class = "Entrez2Ensembl",
    contains = "DataFrame"
)
setValidity(
    Class = "Entrez2Ensembl",
    method = function(object) {
        ok <- validate(
            identical(ncol(object), 2L),
            hasColnames(object),
            hasRows(object),
            all(complete.cases(object))
        )
        if (!isTRUE(ok)) {
            return(ok)
        }
        cols <- c(
            "entrez" = "entrezId",
            "ensembl" = "ensemblId"
        )
        if (!identical(cols, colnames(object))) {
            colnames(object) <- unname(cols)
        }
        ok <- validate(
            is.integer(object[[cols[["entrez"]]]])
        )
        if (!isTRUE(ok)) {
            return(ok)
        }
        TRUE
    }
)



## FIXME Inherit the documentation from AcidGenerics.

#' Gene-to-symbol mappings
#'
#' @details
#' Contains a `DataFrame` with `geneId` and `geneName` columns.
#'
#' @section Genome metadata:
#'
#' We recommend slotting `organism`, `genomeBuild`, and `ensemblRelease` into
#' `metadata()`.
#'
#' @export
#' @note Updated 2021-08-10.
#'
#' @return `Gene2Symbol`.
setClass(
    Class = "Gene2Symbol",
    contains = "DataFrame"
)
setValidity(
    Class = "Gene2Symbol",
    method = function(object) {
        ok <- validate(
            identical(ncol(object), 2L),
            hasColnames(object),
            hasRows(object),
            all(complete.cases(object))
        )
        if (!isTRUE(ok)) {
            return(ok)
        }
        cols <- c(
            "gene" = "geneId",
            "symbol" = "geneName"
        )
        if (!identical(cols, colnames(object))) {
            colnames(object) <- unname(cols)
        }
        ok <- validate(
            is.character(object[[cols[["gene"]]]]),
            isAny(object[[cols[["symbol"]]]], c("character", "factor"))
        )
        if (!isTRUE(ok)) {
            return(ok)
        }
        TRUE
    }
)



## FIXME Inherit the documentation from AcidGenerics.

#' HGNC-to-Ensembl gene identifier mappings
#'
#' @details
#' Contains a `DataFrame` with `hgnc` and `ensembl` columns.
#'
#' @export
#' @note Updated 2021-08-10.
#'
#' @return `HGNC2Ensembl`.
setClass(
    Class = "HGNC2Ensembl",
    contains = "DataFrame"
)
setValidity(
    Class = "HGNC2Ensembl",
    method = function(object) {
        ok <- validate(
            identical(ncol(object), 2L),
            hasColnames(object),
            hasRows(object),
            all(complete.cases(object))
        )
        if (!isTRUE(ok)) {
            return(ok)
        }
        cols <- c(
            "hgnc" = "hgncId",
            "ensembl" = "ensemblId"
        )
        if (!identical(cols, colnames(object))) {
            colnames(object) <- unname(cols)
        }
        ok <- validate(
            is.integer(object[[cols[["hgnc"]]]]),
            hasNoDuplicates(object[[cols[["hgnc"]]]])
        )
        if (!isTRUE(ok)) {
            return(ok)
        }
        TRUE
    }
)



## FIXME Inherit the documentation from AcidGenerics.

#' MGI-to-Ensembl gene identifier mappings
#'
#' @details
#' Contains a `DataFrame` with `mgi` and `ensembl` columns.
#'
#' @export
#' @note Updated 2021-08-10.
#'
#' @return `MGI2Ensembl`.
setClass(
    Class = "MGI2Ensembl",
    contains = "DataFrame"
)
setValidity(
    Class = "MGI2Ensembl",
    method = function(object) {
        ok <- validate(
            identical(ncol(object), 2L),
            hasColnames(object),
            hasRows(object),
            all(complete.cases(object))
        )
        if (!isTRUE(ok)) {
            return(ok)
        }
        cols <- c(
            "mgi" = "mgiId",
            "ensembl" = "ensemblId"
        )
        if (!identical(cols, colnames(object))) {
            colnames(object) <- unname(cols)
        }
        ok <- validate(
            is.integer(object[[cols[["mgi"]]]]),
            hasNoDuplicates(object[[cols[["mgi"]]]])
        )
        if (!isTRUE(ok)) {
            return(ok)
        }
        TRUE
    }
)



## FIXME Inherit the documentation from AcidGenerics.

#' Protein-to-gene mappings
#'
#' @details
#' Contains a `DataFrame` with `proteinId`, `geneId`, and `geneName` columns.
#'
#' @section Genome metadata:
#'
#' We recommend slotting `organism`, `genomeBuild`, and `ensemblRelease` into
#' `metadata()`.
#'
#' @export
#' @note Updated 2021-08-10.
#'
#' @return `Protein2Gene`.
setClass(
    Class = "Protein2Gene",
    contains = "DataFrame"
)
setValidity(
    Class = "Protein2Gene",
    method = function(object) {
        ok <- validate(
            identical(ncol(object), 3L),
            hasColnames(object),
            hasRows(object),
            all(complete.cases(object))
        )
        if (!isTRUE(ok)) {
            return(ok)
        }
        cols <- c("proteinId", "geneId", "geneName")
        if (!identical(cols, colnames(object))) {
            colnames(object) <- camelCase(colnames(object), strict = TRUE)
        }
        ok <- validate(
            identical(cols, colnames(object)),
            all(bapply(X = object, FUN = is.character))
        )
        if (!isTRUE(ok)) {
            return(ok)
        }
        TRUE
    }
)



## FIXME Inherit the documentation from AcidGenerics.

#' Transcript-to-gene identifier mappings
#'
#' @details
#' Contains a `DataFrame` with `txId` and `geneId` columns.
#'
#' @section Genome metadata:
#'
#' We recommend slotting `organism`, `genomeBuild`, and `release` into
#' `metadata`.
#'
#' Ensembl examples:
#' - `organism`: "Homo sapiens".
#' - `genomeBuild`: "GRCh38".
#' - `release`: 100L.
#'
#' @export
#' @note Updated 2021-08-10.
#'
#' @return `Tx2Gene`.
setClass(
    Class = "Tx2Gene",
    contains = "DataFrame"
)
setValidity(
    Class = "Tx2Gene",
    method = function(object) {
        ok <- validate(
            identical(ncol(object), 2L),
            hasColnames(object),
            hasRows(object),
            all(complete.cases(object))
        )
        if (!isTRUE(ok)) {
            return(ok)
        }
        cols <- c(
            "tx" = "txId",
            "gene" = "geneId"
        )
        if (!identical(cols, colnames(object))) {
            colnames(object) <- unname(cols)
        }
        ok <- validate(
            all(vapply(
                X = object,
                FUN = is.character,
                FUN.VALUE = logical(1L)
            )),
            hasNoDuplicates(object[[cols[["tx"]]]])
        )
        if (!isTRUE(ok)) {
            return(ok)
        }
        ok <- validate(
            !any(apply(
                X = object,
                MARGIN = 1L,
                FUN = function(x) {
                    identical(x[[cols[["tx"]]]], x[[cols[["gene"]]]])
                }
            )),
            msg = "Some transcript and gene identifiers are identical."
        )
        if (!isTRUE(ok)) {
            return(ok)
        }
        TRUE
    }
)
