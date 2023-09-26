## Internal validity methods ===================================================

#' Shared Ensembl validity checks
#'
#' @note Updated 2023-04-13.
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



#' Shared GRanges validity checks
#'
#' @note Updated 2023-04-26.
#' @noRd
#'
#' @details
#' Genome build and organism are not defined in minimal FlyBase GTF example.
.grangesValidity <- function(object) {
    ok <- validateClasses(
        object = metadata(object),
        expected = list(
            "date" = "Date",
            "ignoreVersion" = "logical",
            "level" = "character",
            "packageVersion" = "package_version",
            "provider" = "character"
        ),
        subset = TRUE
    )
    if (!isTRUE(ok)) {
        return(ok)
    }
    ok <- validate(
        !isSubset("entrezId", colnames(mcols(object))),
        msg = sprintf(
            "Object contains {.var %s} instead of {.var %s} in {.var %s}.",
            "entrezId", "ncbiGeneId", "mcols"
        )
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



## Genome annotation classes ===================================================

## FIXME Check metadata for date, packageVersion.

#' Ensembl gene annotations
#'
#' @details
#' Contains `GRanges` with Ensembl gene-level annotations.
#'
#' @export
#' @note Updated 2023-04-13.
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



## FIXME Check metadata for date, packageVersion.

#' Ensembl transcript annotations
#'
#' @details
#' Contains `GRanges` with Ensembl transcript-level annotations.
#'
#' @export
#' @note Updated 2023-04-13.
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



## FIXME Check metadata for date, packageVersion.

#' FlyBase gene annotations
#'
#' @details
#' Contains `GRanges` with FlyBase gene-level annotations.
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



## FIXME Check metadata for date, packageVersion.

#' FlyBase transcript annotations
#'
#' @details
#' Contains `GRanges` with FlyBase transcript-level annotations.
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



## FIXME Check metadata for date, packageVersion.

#' GENCODE gene annotations
#'
#' @details
#' Contains `GRanges` with GENCODE gene-level annotations.
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



## FIXME Check metadata for date, packageVersion.

#' GENCODE transcript annotations
#'
#' @details
#' Contains `GRanges` with GENCODE transcript-level annotations.
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



## FIXME Check metadata for date, packageVersion.

#' RefSeq gene annotations
#'
#' @details
#' Contains a `CompressedGRangesList` with RefSeq gene-level annotations.
#'
#' @export
#' @note Updated 2023-04-26.
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



## FIXME Check metadata for date, packageVersion.

#' RefSeq transcript annotations
#'
#' @details
#' Contains a `CompressedGRangesList` with RefSeq transcript-level annotations.
#'
#' @export
#' @note Updated 2023-04-26.
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



## FIXME Check metadata for date, packageVersion.

#' UCSC gene annotations
#'
#' @details
#' Contains `GRanges` with UCSC gene-level annotations.
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



## FIXME Check metadata for date, packageVersion.

#' UCSC transcript annotations
#'
#' @details
#' Contains `GRanges` with UCSC transcript-level annotations.
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



## FIXME Check metadata for date, packageVersion.

#' WormBase gene annotations
#'
#' @details
#' Contains `GRanges` with WormBase gene-level annotations.
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



## FIXME Check metadata for date, packageVersion.

#' WormBase transcript annotations
#'
#' @details
#' Contains `GRanges` with WormBase transcript-level annotations.
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

## FIXME Check metadata for date, packageVersion.

#' Human Genome Organization (HUGO) Gene Nomenclature Committee (HGNC) metadata
#'
#' @export
#' @note Updated 2023-03-01.
#'
#' @return `HGNC`.
setClass(
    Class = "HGNC",
    contains = "DFrame"
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
        ## FIXME Use validate classes here instead.
        cols <- c("hgncId", "ensemblGeneId", "ncbiGeneId")
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



## FIXME Check metadata for date, packageVersion.

#' Mouse Genomic Informatics (MGI) metadata
#'
#' @export
#' @note Updated 2023-09-19.
#'
#' @return `MGI`.
setClass(
    Class = "MGI",
    contains = "DFrame"
)
setValidity(
    Class = "MGI",
    method = function(object) {
        ok <- validate(
            hasColnames(object),
            hasRows(object)
        )
        if (!isTRUE(ok)) {
            return(ok)
        }
        ## FIXME Use validate classes here instead.
        cols <- c(
            "mgiAccessionId",
            "markerType",
            "markerSymbol",
            "markerName",
            "genomeBuild",
            "ncbiGeneId",
            "ncbiGeneChromosome",
            "ncbiGeneStart",
            "ncbiGeneEnd",
            "ncbiGeneStrand",
            "ensemblGeneId",
            "ensemblGeneChromosome",
            "ensemblGeneStart",
            "ensemblGeneEnd",
            "ensemblGeneStrand"
        )
        #ok <- validate(
        #    isSubset(cols, colnames(object)),
        #    is.integer(object[[cols[[1L]]]]),
        #    hasNoDuplicates(object[[cols[[1L]]]])
        #)
        #if (!isTRUE(ok)) {
        #    return(ok)
        #}
        TRUE
    }
)



## FIXME Check metadata for date, packageVersion.

#' NCBI gene history
#'
#' @export
#' @note Updated 2023-09-15.
#'
#' @return `NcbiGeneHistory`.
setClass(
    Class = "NcbiGeneHistory",
    contains = "DFrame"
)
setValidity(
    Class = "NcbiGeneHistory",
    method = function(object) {
        ok <- validate(
            hasRownames(object),
            hasColnames(object)
        )
        if (!isTRUE(ok)) {
            return(ok)
        }
        ok <- validateClasses(
            object = object,
            expected = list(
                "discontinuedGeneId" = "integer",
                "discontinuedSymbol" = "character",
                "discontinueDate" = "Date",
                "geneId" = "integer"
            )
        )
        if (!isTRUE(ok)) {
            return(ok)
        }
        TRUE
    }
)



## FIXME Check metadata for date, packageVersion.

#' NCBI gene identifier information
#'
#' @export
#' @note Updated 2023-03-01.
#'
#' @return `NcbiGeneInfo`.
setClass(
    Class = "NcbiGeneInfo",
    contains = "DFrame"
)
setValidity(
    Class = "NcbiGeneInfo",
    method = function(object) {
        ok <- validate(
            hasColnames(object),
            hasRows(object)
        )
        if (!isTRUE(ok)) {
            return(ok)
        }
        ## FIXME Use validate classes here instead.
        cols <- c(
            "chromosome",
            "dbXrefs",
            "description",
            "featureType",
            "geneId",
            "geneName",
            "geneSynonyms",
            "mapLocation",
            "modificationDate",
            "nomenclatureStatus",
            "otherDesignations",
            "taxonomyId",
            "typeOfGene"
        )
        ok <- validate(
            isSubset(cols, colnames(object))
        )
        if (!isTRUE(ok)) {
            return(ok)
        }
        TRUE
    }
)



## Identifier mapping classes ==================================================

#' @inherit AcidGenerics::EnsemblToNcbi description return title
#' @export
#' @note Updated 2023-09-16.
#'
#' @details
#' Contains a `DFrame` with `ensemblGeneId` and `ncbiGeneId` columns.
setClass(
    Class = "EnsemblToNcbi",
    contains = "DFrame"
)
setValidity(
    Class = "EnsemblToNcbi",
    method = function(object) {
        ok <- validate(
            identical(ncol(object), 2L),
            hasRows(object),
            hasColnames(object),
            all(complete.cases(object))
        )
        if (!isTRUE(ok)) {
            return(ok)
        }
        ok <- validateClasses(
            object = object,
            expected = list(
                "ensemblGeneId" = "character",
                "ncbiGeneId" = "integer"
            )
        )
        if (!isTRUE(ok)) {
            return(ok)
        }
        ok <- validateClasses(
            object = metadata(object),
            expected = list(
                "date" = "Date",
                "format" = "character",
                "packageVersion" = "package_version"
            )
        )
        if (!isTRUE(ok)) {
            return(ok)
        }
        TRUE
    }
)



#' @inherit AcidGenerics::GeneToSymbol description return title
#' @export
#' @note Updated 2023-09-16.
#'
#' @details
#' Contains a `DFrame` with `geneId` and `geneName` columns.
#'
#' @section Genome metadata:
#'
#' We recommend slotting `organism`, `genomeBuild`, and `ensemblRelease` into
#' `metadata()`.
setClass(
    Class = "GeneToSymbol",
    contains = "DFrame"
)
setValidity(
    Class = "GeneToSymbol",
    method = function(object) {
        ok <- validate(
            identical(ncol(object), 2L),
            hasRows(object),
            hasColnames(object),
            all(complete.cases(object))
        )
        if (!isTRUE(ok)) {
            return(ok)
        }
        ok <- validateClasses(
            object = object,
            expected = list(
                "geneId" = "character",
                "geneName" = "character"
            )
        )
        if (!isTRUE(ok)) {
            return(ok)
        }
        ok <- validateClasses(
            object = metadata(object),
            expected = list(
                "date" = "Date",
                "format" = "character",
                "packageVersion" = "package_version"
            ),
            subset = TRUE
        )
        if (!isTRUE(ok)) {
            return(ok)
        }
        TRUE
    }
)



#' Human-to-mouse gene mappings
#'
#' @export
#' @note Updated 2023-09-26.
#'
#' @details
#' Contains a `DFrame` with `"dbClassKey"`, `"humanGeneName"`, `"humanHgncId"`,
#' `"humanNcbiGeneId"`, `"humanOmimGeneId"`, `"mouseGeneName"`, `"mouseMgiId"`,
#' `"mouseNcbiGeneId"` columns.
#'
#' @return `HumanToMouse`.
setClass(
    Class = "HumanToMouse",
    contains = "DFrame"
)
setValidity(
    Class = "HumanToMouse",
    method = function(object) {
        ok <- validate(
            hasRows(object),
            hasColnames(object)
        )
        if (!isTRUE(ok)) {
            return(ok)
        }
        ok <- validateClasses(
            object = object,
            expected = list(
                "dbClassKey" = "integer",
                "humanGeneName" = "character",
                "humanHgncId" = "integer",
                "humanNcbiGeneId" = "integer",
                "humanOmimGeneId" = "integer",
                "mouseGeneName" = "character",
                "mouseMgiId" = "integer",
                "mouseNcbiGeneId" = "integer"
            )
        )
        if (!isTRUE(ok)) {
            return(ok)
        }
        ok <- validateClasses(
            object = metadata(object),
            expected = list(
                "date" = "Date",
                "packageVersion" = "package_version",
                "unique" = "logical"
            ),
            subset = TRUE
        )
        if (!isTRUE(ok)) {
            return(ok)
        }
        TRUE
    }
)



#' @inherit AcidGenerics::NcbiToEnsembl description return title
#' @export
#' @note Updated 2023-09-16.
#'
#' @details
#' Contains a `DFrame` with `ncbiGeneId` and `ensemblGeneId` columns.
setClass(
    Class = "NcbiToEnsembl",
    contains = "DFrame"
)
setValidity(
    Class = "NcbiToEnsembl",
    method = function(object) {
        ok <- validate(
            identical(ncol(object), 2L),
            hasRows(object),
            hasColnames(object),
            all(complete.cases(object))
        )
        if (!isTRUE(ok)) {
            return(ok)
        }
        ok <- validateClasses(
            object = object,
            expected = list(
                "ncbiGeneId" = "integer",
                "ensemblGeneId" = "character"
            )
        )
        if (!isTRUE(ok)) {
            return(ok)
        }
        ok <- validateClasses(
            object = metadata(object),
            expected = list(
                "date" = "Date",
                "format" = "character",
                "packageVersion" = "package_version"
            )
        )
        if (!isTRUE(ok)) {
            return(ok)
        }
        TRUE
    }
)



#' Protein-to-gene mappings
#'
#' @details
#' Contains a `DFrame` with `proteinId`, `geneId`, and `geneName` columns.
#'
#' @section Genome metadata:
#'
#' We recommend slotting `organism`, `genomeBuild`, and `ensemblRelease` into
#' `metadata()`.
#'
#' @export
#' @note Updated 2023-09-16.
#'
#' @return `ProteinToGene`.
setClass(
    Class = "ProteinToGene",
    contains = "DFrame"
)
setValidity(
    Class = "ProteinToGene",
    method = function(object) {
        ok <- validate(
            identical(ncol(object), 3L),
            hasRows(object),
            hasColnames(object),
            all(complete.cases(object))
        )
        if (!isTRUE(ok)) {
            return(ok)
        }
        ok <- validateClasses(
            object = object,
            expected = list(
                "proteinId" = "character",
                "geneId" = "character",
                "geneName" = "character"
            )
        )
        if (!isTRUE(ok)) {
            return(ok)
        }
        ok <- validateClasses(
            object = metadata(object),
            expected = list(
                "date" = "Date",
                "packageVersion" = "package_version"
            ),
            subset = TRUE
        )
        if (!isTRUE(ok)) {
            return(ok)
        }
        TRUE
    }
)



## FIXME Check metadata for date, packageVersion.

#' @inherit AcidGenerics::TxToGene description return title
#' @export
#' @note Updated 2023-09-16.
#'
#' @details
#' Contains a `DFrame` with `txId` and `geneId` columns.
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
setClass(
    Class = "TxToGene",
    contains = "DFrame"
)
setValidity(
    Class = "TxToGene",
    method = function(object) {
        ok <- validate(
            identical(ncol(object), 2L),
            hasRows(object),
            hasColnames(object),
            all(complete.cases(object))
        )
        if (!isTRUE(ok)) {
            return(ok)
        }
        ## FIXME Use validate classes here instead.
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
        ## This check for identical transcript and gene identifiers is not
        ## compatible with the sacCer3 genome, unfortunately. Disabling this
        ## check until we can come up with an improved approach.
        ## See related issue:
        ## https://github.com/bcbio/bcbio-nextgen/issues/3565
        ## > ok <- validate(
        ## >     !any(apply(
        ## >         X = object,
        ## >         MARGIN = 1L,
        ## >         FUN = function(x) {
        ## >             identical(x[[cols[["tx"]]]], x[[cols[["gene"]]]])
        ## >         }
        ## >     )),
        ## >     msg = "Some transcript and gene identifiers are identical."
        ## > )
        ## > if (!isTRUE(ok)) return(ok)
        ##
        ## FIXME Add metadata check here.
        TRUE
    }
)
