## FIXME Need to harden validity checks for classed GRanges.



## Genome annotation classes ===================================================
#' Shared Ensembl validity checks
#'
#' @note Updated 2021-01-22.
#' @noRd
.ensemblValidity <- function(object) {
    validate(
        ## GRCh37 currently will fail these checks.
        ## > !all(is.na(seqlengths(object))),
        ## > !all(is.na(seqnames(object))),
        ## > !all(is.na(genome(object))),
        identical(
            x = colnames(mcols(object)),
            y = camelCase(colnames(mcols(object)), strict = TRUE)
        ),
        isString(metadata(object)[["genomeBuild"]]),
        isFlag(metadata(object)[["ignoreVersion"]]),
        isOrganism(metadata(object)[["organism"]])
        ## Not necessarily defined in GFF file.
        ## > isInt(metadata(object)[["release"]])
    )
}



#' Shared FlyBase validity checks
#'
#' @note Updated 2021-01-25.
#' @noRd
.flybaseValidity <- function(object) {
    validate(
        identical(metadata(object)[["provider"]], "FlyBase")
    )
}



#' Shared GENCODE validity checks
#'
#' @note Updated 2021-01-25.
#' @noRd
.gencodeValidity <- function(object) {
    validate(
        identical(metadata(object)[["provider"]], "GENCODE")
    )
}



#' Shared RefSeq validity checks
#'
#' @note Updated 2021-01-25.
#' @noRd
.refseqValidity <- function(object) {
    validate(
        identical(metadata(object)[["provider"]], "RefSeq")
    )
}



#' Shared UCSC validity checks
#'
#' @note Updated 2021-01-25.
#' @noRd
.ucscValidity <- function(object) {
    validate(
        identical(metadata(object)[["provider"]], "UCSC")
    )
}



#' Shared WormBase validity checks
#'
#' @note Updated 2021-01-25.
#' @noRd
.wormbaseValidity <- function(object) {
    validate(
        identical(metadata(object)[["provider"]], "WormBase")
    )
}



#' Ensembl gene annotations
#'
#' @details
#' Contains a `GRanges` with Ensembl gene-level annotations.
#'
#' @export
#' @note Updated 2021-01-22.
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
        if (!isTRUE(ok)) return(ok)
        ok <- validate(
            any(grepl(
                pattern = "^ENS([A-Z]+)?G[0-9]{11}(\\.[0-9]+)?$",
                x = names(object)
            )),
            identical(metadata(object)[["level"]], "genes")
        )
        if (!isTRUE(ok)) return(ok)
        TRUE
    }
)



#' Ensembl transcript annotations
#'
#' @details
#' Contains a `GRanges` with Ensembl transcript-level annotations.
#'
#' @export
#' @note Updated 2021-01-22.
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
        if (!isTRUE(ok)) return(ok)
        ok <- validate(
            any(grepl(
                pattern = "^ENS([A-Z]+)?T[0-9]{11}(\\.[0-9]+)?$",
                x = names(object)
            )),
            identical(metadata(object)[["level"]], "transcripts")
        )
        if (!isTRUE(ok)) return(ok)
        TRUE
    }
)



#' FlyBase gene annotations
#'
#' @details
#' Contains a `GRanges` with FlyBase gene-level annotations.
#'
#' @export
#' @note Updated 2021-01-25.
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
        if (!isTRUE(ok)) return(ok)
        TRUE
    }
)



#' FlyBase transcript annotations
#'
#' @details
#' Contains a `GRanges` with FlyBase transcript-level annotations.
#'
#' @export
#' @note Updated 2021-01-25.
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
        if (!isTRUE(ok)) return(ok)
        TRUE
    }
)



#' GENCODE gene annotations
#'
#' @details
#' Contains a `GRanges` with GENCODE gene-level annotations.
#'
#' @export
#' @note Updated 2021-01-25.
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
        if (!isTRUE(ok)) return(ok)
        TRUE
    }
)



#' GENCODE transcript annotations
#'
#' @details
#' Contains a `GRanges` with GENCODE transcript-level annotations.
#'
#' @export
#' @note Updated 2021-01-25.
#'
#' @return `GencodeTranscripts`.
setClass(
    Class = "GencodeTranscripts",
    contains = "EnsemblTranscripts"
)
setValidity(
    Class = "GencodeTranscripts",
    method = function(object) {
        ok <- .gencodeValidity(object)
        if (!isTRUE(ok)) return(ok)
        TRUE
    }
)



#' RefSeq gene annotations
#'
#' @details
#' Contains a `GRanges` with RefSeq gene-level annotations.
#'
#' @export
#' @note Updated 2021-01-26.
#'
#' @return `RefSeqGenes`.
setClass(
    Class = "RefSeqGenes",
    contains = "GRanges"
)
setValidity(
    Class = "CompressedGRangesList",
    method = function(object) {
        ## FIXME Is this not setting metadata correctly???
        ## > ok <- .refseqValidity(object)
        ## > if (!isTRUE(ok)) return(ok)
        TRUE
    }
)



#' RefSeq transcript annotations
#'
#' @details
#' Contains a `GRanges` with RefSeq transcript-level annotations.
#'
#' @export
#' @note Updated 2021-01-26.
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
        if (!isTRUE(ok)) return(ok)
        TRUE
    }
)



#' UCSC gene annotations
#'
#' @details
#' Contains a `GRanges` with UCSC gene-level annotations.
#'
#' @export
#' @note Updated 2021-01-25.
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
        if (!isTRUE(ok)) return(ok)
        TRUE
    }
)



#' UCSC transcript annotations
#'
#' @details
#' Contains a `GRanges` with UCSC transcript-level annotations.
#'
#' @export
#' @note Updated 2021-01-25.
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
        if (!isTRUE(ok)) return(ok)
        TRUE
    }
)



#' WormBase gene annotations
#'
#' @details
#' Contains a `GRanges` with WormBase gene-level annotations.
#'
#' @export
#' @note Updated 2021-01-25.
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
        if (!isTRUE(ok)) return(ok)
        TRUE
    }
)



#' WormBase transcript annotations
#'
#' @details
#' Contains a `GRanges` with WormBase transcript-level annotations.
#'
#' @export
#' @note Updated 2021-01-25.
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
        if (!isTRUE(ok)) return(ok)
        TRUE
    }
)



## Identifier classes ==========================================================
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
            isSubset(c("hgncId", "ensemblGeneId"), colnames(object)),
            is.integer(object[["hgncId"]]),
            hasNoDuplicates(object[["hgncId"]])
        )
    }
)



## Identifier mapping classes ==================================================
#' Ensembl-to-Entrez gene identifier mappings
#'
#' @details
#' Contains a `DataFrame` with `ensemblId` and `entrezId` columns.
#'
#' @export
#' @note Updated 2021-01-18.
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
            identical(c("ensemblId", "entrezId"), colnames(object)),
            is.integer(object[["entrezId"]])
        )
    }
)



#' Entrez-to-Ensembl gene identifier mappings
#'
#' @inherit Ensembl2Entrez-class details
#'
#' @export
#' @note Updated 2021-01-18.
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
            identical(c("entrezId", "ensemblId"), colnames(object)),
            is.integer(object[["entrezId"]])
        )
    }
)



## FIXME THIS NEEDS TO SLOT `ignoreVersions` in metadata.

#' Gene-to-symbol mappings
#'
#' @details
#' Contains a `DataFrame` with `geneId` and `geneName` columns.
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
            identical(colnames(object), c("geneId", "geneName")),
            is.character(object[["geneId"]]),
            isAny(object[["geneName"]], c("character", "factor"))
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
#' Contains a `DataFrame` with `proteinId`, `geneId`, and `geneName` columns.
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
                y = c("proteinId", "geneId", "geneName")
            ),
            all(bapply(X = object, FUN = is.character))
        )
    }
)



## FIXME THIS NEEDS TO SLOT `ignoreVersions` in metadata.

#' Transcript-to-gene identifier mappings
#'
#' @details
#' Contains a `DataFrame` with `txId` and `geneId` columns.
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
            identical(colnames(object), c("txId", "geneId")),
            all(vapply(
                X = object,
                FUN = is.character,
                FUN.VALUE = logical(1L)
            )),
            hasNoDuplicates(object[["txId"]])
        )
    }
)
