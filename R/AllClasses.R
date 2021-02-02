## Genome annotation classes ===================================================
#' Shared GRanges validity checks
#'
#' @note Note that genome build and organism are not defined in minimal FlyBase
#'   GTF example.
#' @note Updated 2021-02-01.
#' @noRd
.grangesValidity <- function(object) {
    if (is(object, "GRangesList")) {
        gr <- object[[1L]]
    } else {
        gr <- object
    }
    ok <- validate(
        identical(
            x = colnames(mcols(gr)),
            y = camelCase(colnames(mcols(gr)), strict = TRUE)
        )
    )
    if (!isTRUE(ok)) return(ok)
    ok <- validateClasses(
        object = metadata(object),
        expected = list(
            "acidGenomes" = "package_version",
            "date" = "Date",
            "ignoreVersion" = "logical",
            "level" = "character",
            "provider" = "character",
            "synonyms" = "logical"
        ),
        subset = TRUE
    )
    if (!isTRUE(ok)) return(ok)
    TRUE
}



#' Shared Ensembl validity checks
#'
#' @note Updated 2021-02-01.
#' @noRd
.ensemblValidity <- function(object) {
    ok <- .grangesValidity(object)
    if (!isTRUE(ok)) return(ok)
    ok <- validate(
        identical(metadata(object)[["provider"]], "Ensembl")
    )
    if (!isTRUE(ok)) return(ok)
    ok <- validateClasses(
        object = metadata(object),
        expected = list(
            "genomeBuild" = "character",
            "organism" = "character"
        ),
        subset = TRUE
    )
    if (!isTRUE(ok)) return(ok)
    ## Can't always get the release version from GFF file, so just check
    ## for ensembldb return.
    if (isSubset("ensembldb", names(metadata(object)))) {
        ok <- validate(
            is.integer(metadata(object)[["release"]])
        )
        if (!isTRUE(ok)) return(ok)
    }
    TRUE
}



#' Shared FlyBase validity checks
#'
#' @note Updated 2021-01-30.
#' @noRd
.flybaseValidity <- function(object) {
    ok <- .grangesValidity(object)
    if (!isTRUE(ok)) return(ok)
    ok <- validate(
        identical(metadata(object)[["provider"]], "FlyBase")
    )
    if (!isTRUE(ok)) return(ok)
    TRUE
}



#' Shared GENCODE validity checks
#'
#' @note Updated 2021-02-01.
#' @noRd
.gencodeValidity <- function(object) {
    ok <- .grangesValidity(object)
    if (!isTRUE(ok)) return(ok)
    ok <- validate(
        identical(metadata(object)[["provider"]], "GENCODE")
    )
    if (!isTRUE(ok)) return(ok)
    ok <- validateClasses(
        object = metadata(object),
        expected = list(
            "genomeBuild" = "character",
            "organism" = "character"
        ),
        subset = TRUE
    )
    if (!isTRUE(ok)) return(ok)
    TRUE
}



#' Shared RefSeq validity checks
#'
#' @note Updated 2021-02-01.
#' @noRd
.refseqValidity <- function(object) {
    ok <- .grangesValidity(object)
    if (!isTRUE(ok)) return(ok)
    ok <- validate(
        identical(metadata(object)[["provider"]], "RefSeq")
    )
    if (!isTRUE(ok)) return(ok)
    ok <- validateClasses(
        object = metadata(object),
        expected = list(
            "genomeBuild" = "character",
            "organism" = "character"
        ),
        subset = TRUE
    )
    if (!isTRUE(ok)) return(ok)
    TRUE
}



#' Shared UCSC validity checks
#'
#' @note Updated 2021-02-01.
#' @noRd
.ucscValidity <- function(object) {
    ok <- .grangesValidity(object)
    if (!isTRUE(ok)) return(ok)
    ok <- validate(
        identical(metadata(object)[["provider"]], "UCSC")
    )
    if (!isTRUE(ok)) return(ok)
    ok <- validateClasses(
        object = metadata(object),
        expected = list(
            "genomeBuild" = "character",
            "organism" = "character"
        ),
        subset = TRUE
    )
    if (!isTRUE(ok)) return(ok)
    TRUE
}



#' Shared WormBase validity checks
#'
#' @note Updated 2021-02-01.
#' @noRd
.wormbaseValidity <- function(object) {
    ok <- .grangesValidity(object)
    if (!isTRUE(ok)) return(ok)
    ok <- validate(
        identical(metadata(object)[["provider"]], "WormBase")
    )
    if (!isTRUE(ok)) return(ok)
    TRUE
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
            allAreMatchingRegex(
                pattern = paste0(
                    "^",
                    "(",
                    "ENS([A-Z]+)?G[0-9]{11}",
                    "|",
                    "LRG_[0-9]+",
                    ")",
                    "(\\.[0-9]+)?",
                    "$"
                ),
                x = names(object)
            ),
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
            allAreMatchingRegex(
                pattern = paste0(
                    "^",
                    "(",
                    "ENS([A-Z]+)?T[0-9]{11}",
                    "|",
                    "LRG_[0-9]+t[0-9]+(-[0-9]+)?",
                    ")",
                    "(\\.[0-9]+)?",
                    "$"
                ),
                x = names(object)
            ),
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
    contains = "GRanges"
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
#' Contains a `GRangesList` with RefSeq gene-level annotations.
#'
#' @export
#' @note Updated 2021-01-30.
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
        if (!isTRUE(ok)) return(ok)
        TRUE
    }
)



#' RefSeq transcript annotations
#'
#' @details
#' Contains a `GRangesList` with RefSeq transcript-level annotations.
#'
#' @export
#' @note Updated 2021-01-30.
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
#' @note Updated 2021-02-02.
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
        if (!isTRUE(ok)) return(ok)
        cols <- c("hgncId", "ensemblGeneId")
        if (!isSubset(cols, colnames(object))) {
            colnames(object) <- camelCase(colnames(object), strict = TRUE)
        }
        ok <- validate(
            isSubset(cols, colnames(object)),
            is.integer(object[[cols[[1L]]]]),
            hasNoDuplicates(object[[cols[[1L]]]])
        )
        if (!isTRUE(ok)) return(ok)
        TRUE
    }
)



## Identifier mapping classes ==================================================
#' Ensembl-to-Entrez gene identifier mappings
#'
#' @details
#' Contains a `DataFrame` with `ensemblId` and `entrezId` columns.
#'
#' @export
#' @note Updated 2021-02-02.
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
            hasRows(object)
        )
        if (!isTRUE(ok)) return(ok)
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
        if (!isTRUE(ok)) return(ok)
        TRUE
    }
)



#' Entrez-to-Ensembl gene identifier mappings
#'
#' @inherit Ensembl2Entrez-class details
#'
#' @export
#' @note Updated 2021-02-02.
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
            hasRows(object)
        )
        if (!isTRUE(ok)) return(ok)
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
        if (!isTRUE(ok)) return(ok)
        TRUE
    }
)



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
#' @note Updated 2021-02-02.
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
            hasRows(object)
        )
        if (!isTRUE(ok)) return(ok)
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
        if (!isTRUE(ok)) return(ok)
        TRUE
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
#' @note Updated 2021-02-01.
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
            identical(colnames(object), c("mgiId", "ensemblId")),
            is.integer(object[["mgiId"]]),
            hasNoDuplicates(object[["mgiId"]])
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
#' @note Updated 2021-01-29.
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
            hasRows(object),
            identical(ncol(object), 2L)
        )
        if (!isTRUE(ok)) return(ok)
        ## Note that "transcriptId" is allowed for legacy compatibility.
        ok <- validate(
            isSubset(
                x = camelCase(colnames(object)[[1L]], strict = TRUE),
                y = c("transcriptId", "txId")
            ),
            identical(
                x = camelCase(colnames(object)[[2L]], strict = TRUE),
                y = "geneId"
            ),
            msg = "Column names are invalid. Use 'txId' and 'geneId'."
        )
        if (!isTRUE(ok)) return(ok)
        ok <- validate(
            all(vapply(
                X = object,
                FUN = is.character,
                FUN.VALUE = logical(1L)
            )),
            hasNoDuplicates(object[[1L]])
        )
        if (!isTRUE(ok)) return(ok)
        ok <- validate(
            !any(apply(
                X = object,
                MARGIN = 1L,
                FUN = function(x) {
                    identical(x[[1L]], x[[2L]])
                }
            )),
            msg = "Some transcript and gene identifiers are identical."
        )
        if (!isTRUE(ok)) return(ok)
        TRUE
    }
)
