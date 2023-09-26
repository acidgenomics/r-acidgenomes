## Internal validity methods ===================================================

#' Shared Ensembl validity checks
#'
#' @note Updated 2023-09-26.
#' @noRd
.validateEnsembl <- function(object) {
    ok <- .validateGRanges(object)
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
        expected = list("organism" = "character"),
        subset = TRUE
    )
    if (!isTRUE(ok)) {
        return(ok)
    }
    ## Can't always get the release version from GFF file, so just check
    ## for ensembldb return.
    if (isSubset("ensembldb", names(metadata(object)))) {
        ok <- validateClasses(
            object = metadata(object),
            expected = list("release" = "integer"),
            subset = TRUE
        )
        if (!isTRUE(ok)) {
            return(ok)
        }
    }
    TRUE
}



#' Shared FlyBase validity checks
#'
#' @note Updated 2023-09-26.
#' @noRd
.validateFlybase <- function(object) {
    ok <- .validateGRanges(object)
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
#' @note Updated 2023-09-26.
#' @noRd
.validateGencode <- function(object) {
    ok <- .validateGRanges(object)
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
#' @note Updated 2023-09-26.
#' @noRd
#'
#' @details
#' Genome build and organism are not defined in minimal FlyBase GTF example.
.validateGRanges <- function(object) {
    ok <- .validateMetadata(object)
    if (!isTRUE(ok)) {
        return(ok)
    }
    ok <- validateClasses(
        object = metadata(object),
        expected = list(
            "ignoreVersion" = "logical",
            "level" = "character",
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



#' Shared metadata validity checks
#'
#' @note Updated 2023-09-26.
#' @noRd
.validateMetadata <- function(object) {
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



#' Shared RefSeq validity checks
#'
#' @note Updated 2023-09-26.
#' @noRd
.validateRefseq <- function(object) {
    ok <- .validateGRanges(object)
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
#' @note Updated 2023-09-26.
#' @noRd
.validateUcsc <- function(object) {
    ok <- .validateGRanges(object)
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
#' @note Updated 2023-09-26.
#' @noRd
.validateWormbase <- function(object) {
    ok <- .validateGRanges(object)
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

#' Ensembl gene annotations
#'
#' @details
#' Contains `GRanges` with Ensembl gene-level annotations.
#'
#' @export
#' @note Updated 2023-09-26.
#'
#' @return `EnsemblGenes`.
setClass(
    Class = "EnsemblGenes",
    contains = "GRanges"
)
setValidity(
    Class = "EnsemblGenes",
    method = function(object) {
        ok <- .validateEnsembl(object)
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
#' Contains `GRanges` with Ensembl transcript-level annotations.
#'
#' @export
#' @note Updated 2023-09-26.
#'
#' @return `EnsemblTranscripts`.
setClass(
    Class = "EnsemblTranscripts",
    contains = "GRanges"
)
setValidity(
    Class = "EnsemblTranscripts",
    method = function(object) {
        ok <- .validateEnsembl(object)
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
#' Contains `GRanges` with FlyBase gene-level annotations.
#'
#' @export
#' @note Updated 2023-09-26.
#'
#' @return `FlyBaseGenes`.
setClass(
    Class = "FlyBaseGenes",
    contains = "GRanges"
)
setValidity(
    Class = "FlyBaseGenes",
    method = function(object) {
        ok <- .validateFlybase(object)
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



#' FlyBase transcript annotations
#'
#' @details
#' Contains `GRanges` with FlyBase transcript-level annotations.
#'
#' @export
#' @note Updated 2023-09-26.
#'
#' @return `FlyBaseTranscripts`.
setClass(
    Class = "FlyBaseTranscripts",
    contains = "GRanges"
)
setValidity(
    Class = "FlyBaseTranscripts",
    method = function(object) {
        ok <- .validateFlybase(object)
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



#' GENCODE gene annotations
#'
#' @details
#' Contains `GRanges` with GENCODE gene-level annotations.
#'
#' @export
#' @note Updated 2023-09-26.
#'
#' @return `GencodeGenes`.
setClass(
    Class = "GencodeGenes",
    contains = "GRanges"
)
setValidity(
    Class = "GencodeGenes",
    method = function(object) {
        ok <- .validateGencode(object)
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



#' GENCODE transcript annotations
#'
#' @details
#' Contains `GRanges` with GENCODE transcript-level annotations.
#'
#' @export
#' @note Updated 2023-09-26.
#'
#' @return `GencodeTranscripts`.
setClass(
    Class = "GencodeTranscripts",
    contains = "GRanges"
)
setValidity(
    Class = "GencodeTranscripts",
    method = function(object) {
        ok <- .validateGencode(object)
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



#' RefSeq gene annotations
#'
#' @details
#' Contains a `CompressedGRangesList` with RefSeq gene-level annotations.
#'
#' @export
#' @note Updated 2023-09-26.
#'
#' @return `RefSeqGenes`.
setClass(
    Class = "RefSeqGenes",
    contains = "CompressedGRangesList"
)
setValidity(
    Class = "RefSeqGenes",
    method = function(object) {
        ok <- .validateRefseq(object)
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



#' RefSeq transcript annotations
#'
#' @details
#' Contains a `CompressedGRangesList` with RefSeq transcript-level annotations.
#'
#' @export
#' @note Updated 2023-09-26.
#'
#' @return `RefSeqTranscripts`.
setClass(
    Class = "RefSeqTranscripts",
    contains = "CompressedGRangesList"
)
setValidity(
    Class = "RefSeqTranscripts",
    method = function(object) {
        ok <- .validateRefseq(object)
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



#' UCSC gene annotations
#'
#' @details
#' Contains `GRanges` with UCSC gene-level annotations.
#'
#' @export
#' @note Updated 2023-09-26.
#'
#' @return `UCSCGenes`.
setClass(
    Class = "UCSCGenes",
    contains = "GRanges"
)
setValidity(
    Class = "UCSCGenes",
    method = function(object) {
        ok <- .validateUcsc(object)
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



#' UCSC transcript annotations
#'
#' @details
#' Contains `GRanges` with UCSC transcript-level annotations.
#'
#' @export
#' @note Updated 2023-09-26.
#'
#' @return `UCSCTranscripts`.
setClass(
    Class = "UCSCTranscripts",
    contains = "GRanges"
)
setValidity(
    Class = "UCSCTranscripts",
    method = function(object) {
        ok <- .validateUcsc(object)
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



#' WormBase gene annotations
#'
#' @details
#' Contains `GRanges` with WormBase gene-level annotations.
#'
#' @export
#' @note Updated 2023-09-26.
#'
#' @return `WormBaseGenes`.
setClass(
    Class = "WormBaseGenes",
    contains = "GRanges"
)
setValidity(
    Class = "WormBaseGenes",
    method = function(object) {
        ok <- .validateWormbase(object)
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



#' WormBase transcript annotations
#'
#' @details
#' Contains `GRanges` with WormBase transcript-level annotations.
#'
#' @export
#' @note Updated 2023-09-26.
#'
#' @return `WormBaseTranscripts`.
setClass(
    Class = "WormBaseTranscripts",
    contains = "GRanges"
)
setValidity(
    Class = "WormBaseTranscripts",
    method = function(object) {
        ok <- .validateWormbase(object)
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



## Identifier classes ==========================================================

#' Human Genome Organization (HUGO) Gene Nomenclature Committee (HGNC) metadata
#'
#' @export
#' @note Updated 2023-09-26.
#'
#' @return `HGNC`.
setClass(
    Class = "HGNC",
    contains = "DFrame"
)
setValidity(
    Class = "HGNC",
    method = function(object) {
        ok <- validateClasses(
            object = object,
            expected = list(
                "agr" = "character",
                "aliasName" = "CompressedCharacterList",
                "aliasSymbol" = "CompressedCharacterList",
                "bioparadigmsSlc" = "character",
                "ccdsId" = "CompressedCharacterList",
                "cd" = "character",
                "cosmic" = "character",
                "dateApprovedReserved" = "Date",
                "dateModified" = "Date",
                "dateNameChanged" = "Date",
                "dateSymbolChanged" = "Date",
                "ena" = "CompressedCharacterList",
                "ensemblGeneId" = "character",
                "enzymeId" = "CompressedCharacterList",
                "gencc" = "character",
                "geneGroup" = "CompressedCharacterList",
                "geneGroupId" = "CompressedCharacterList",
                "gtrnadb" = "character",
                "hgncId" = "integer",
                "homeodb" = "numeric",
                "hordeId" = "character",
                "imgt" = "character",
                "intermediateFilamentDb" = "logical",
                "iuphar" = "character",
                "kznfGeneCatalog" = "logical",
                "lncipedia" = "character",
                "lncrnadb" = "character",
                "location" = "character",
                "locationSortable" = "character",
                "locusGroup" = "character",
                "locusType" = "character",
                "lsdb" = "CompressedCharacterList",
                "mamitTrnadb" = "numeric",
                "maneSelect" = "CompressedCharacterList",
                "merops" = "character",
                "mgdId" = "CompressedCharacterList",
                "mirbase" = "character",
                "name" = "character",
                "ncbiGeneId" = "integer",
                "omimId" = "CompressedCharacterList",
                "orphanet" = "numeric",
                "prevName" = "CompressedCharacterList",
                "prevSymbol" = "CompressedCharacterList",
                "pseudogeneOrg" = "character",
                "pubmedId" = "CompressedCharacterList",
                "refseqAccession" = "CompressedCharacterList",
                "rgdId" = "CompressedCharacterList",
                "rnaCentralIds" = "logical",
                "snornabase" = "character",
                "status" = "character",
                "symbol" = "character",
                "ucscId" = "character",
                "uniprotIds" = "CompressedCharacterList",
                "vegaId" = "character"
            )
        )
        ok <- validate(hasRows(object))
        if (!isTRUE(ok)) {
            return(ok)
        }
        ok <- .validateMetadata(object)
        if (!isTRUE(ok)) {
            return(ok)
        }
        TRUE
    }
)



#' Mouse Genomic Informatics (MGI) metadata
#'
#' @export
#' @note Updated 2023-09-26.
#'
#' @return `MGI`.
setClass(
    Class = "MGI",
    contains = "DFrame"
)
setValidity(
    Class = "MGI",
    method = function(object) {
        ok <- validateClasses(
            object = object,
            expected = list(
                "ensemblGeneChromosome" = "Rle",
                "ensemblGeneEnd" = "Rle",
                "ensemblGeneId" = "Rle",
                "ensemblGeneStart" = "Rle",
                "ensemblGeneStrand" = "Rle",
                "genomeBuild" = "Rle",
                "markerName" = "Rle",
                "markerSymbol" = "Rle",
                "markerType" = "Rle",
                "mgiAccessionId" = "Rle",
                "ncbiGeneChromosome" = "Rle",
                "ncbiGeneEnd" = "Rle",
                "ncbiGeneId" = "Rle",
                "ncbiGeneStart" = "Rle",
                "ncbiGeneStrand" = "Rle"
            )
        )
        if (!isTRUE(ok)) {
            return(ok)
        }
        ok <- validate(hasRows(object))
        if (!isTRUE(ok)) {
            return(ok)
        }
        ok <- .validateMetadata(object)
        if (!isTRUE(ok)) {
            return(ok)
        }
        ok <- validateClasses(
            object = metadata(object),
            expected = list("url" = "character"),
            subset = TRUE
        )
        if (!isTRUE(ok)) {
            return(ok)
        }
        TRUE
    }
)



#' NCBI gene history
#'
#' @export
#' @note Updated 2023-09-26.
#'
#' @return `NcbiGeneHistory`.
setClass(
    Class = "NcbiGeneHistory",
    contains = "DFrame"
)
setValidity(
    Class = "NcbiGeneHistory",
    method = function(object) {
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
        ok <- validate(hasRownames(object))
        if (!isTRUE(ok)) {
            return(ok)
        }
        ok <- .validateMetadata(object)
        if (!isTRUE(ok)) {
            return(ok)
        }
        TRUE
    }
)



#' NCBI gene identifier information
#'
#' @export
#' @note Updated 2023-09-26.
#'
#' @return `NcbiGeneInfo`.
setClass(
    Class = "NcbiGeneInfo",
    contains = "DFrame"
)
setValidity(
    Class = "NcbiGeneInfo",
    method = function(object) {
        ok <- validateClasses(
            object = object,
            expected = list(
                "chromosome" = "Rle",
                "dbXrefs" = "CompressedCharacterList",
                "description" = "Rle",
                "featureType" = "Rle",
                "geneId" = "Rle",
                "geneName" = "Rle",
                "geneSynonyms" = "CompressedCharacterList",
                "mapLocation" = "Rle",
                "modificationDate" = "Date",
                "nomenclatureStatus" = "Rle",
                "otherDesignations" = "CompressedCharacterList",
                "taxonomyId" = "Rle",
                "typeOfGene" = "Rle"
            )
        )
        if (!isTRUE(ok)) {
            return(ok)
        }
        ok <- .validateMetadata(object)
        if (!isTRUE(ok)) {
            return(ok)
        }
        TRUE
    }
)



## Identifier mapping classes ==================================================

#' @inherit AcidGenerics::EnsemblToNcbi description return title
#' @export
#' @note Updated 2023-09-26.
#'
#' @details
#' Contains a `DFrame` with `"ensemblGeneId"` and `"ncbiGeneId"` columns.
setClass(
    Class = "EnsemblToNcbi",
    contains = "DFrame"
)
setValidity(
    Class = "EnsemblToNcbi",
    method = function(object) {
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
        ok <- validate(
            hasRows(object),
            all(complete.cases(object))
        )
        if (!isTRUE(ok)) {
            return(ok)
        }
        ok <- .validateMetadata(object)
        if (!isTRUE(ok)) {
            return(ok)
        }
        ok <- validateClasses(
            object = metadata(object),
            expected = list("format" = "character"),
            isSubset = TRUE
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
#' Contains a `DFrame` with `"geneId"` and `"geneName"` columns.
#'
#' @section Genome metadata:
#'
#' We recommend slotting `"organism"`, `"genomeBuild"`, and `"ensemblRelease"`
#' into `metadata()`.
setClass(
    Class = "GeneToSymbol",
    contains = "DFrame"
)
setValidity(
    Class = "GeneToSymbol",
    method = function(object) {
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
        ok <- validate(
            hasRows(object),
            all(complete.cases(object))
        )
        if (!isTRUE(ok)) {
            return(ok)
        }
        ok <- .validateMetadata(object)
        if (!isTRUE(ok)) {
            return(ok)
        }
        ok <- validateClasses(
            object = metadata(object),
            expected = list("format" = "character"),
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
        ok <- validate(hasRows(object))
        if (!isTRUE(ok)) {
            return(ok)
        }
        ok <- .validateMetadata(object)
        if (!isTRUE(ok)) {
            return(ok)
        }
        ok <- validateClasses(
            object = metadata(object),
            expected = list(
                "unique" = "logical",
                "url" = "character"
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
#' @note Updated 2023-09-26.
#'
#' @details
#' Contains a `DFrame` with `"ncbiGeneId"` and `"ensemblGeneId"` columns.
setClass(
    Class = "NcbiToEnsembl",
    contains = "DFrame"
)
setValidity(
    Class = "NcbiToEnsembl",
    method = function(object) {
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
        ok <- validate(
            hasRows(object),
            all(complete.cases(object))
        )
        if (!isTRUE(ok)) {
            return(ok)
        }
        ok <- .validateMetadata(object)
        if (!isTRUE(ok)) {
            return(ok)
        }
        ok <- validateClasses(
            object = metadata(object),
            expected = list("format" = "character"),
            subset = TRUE
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
#' Contains a `DFrame` with `"proteinId"`, `"geneId"`, and `"geneName"` columns.
#'
#' @section Genome metadata:
#'
#' We recommend slotting `"organism"`, `"genomeBuild"`, and `"ensemblRelease"`
#' into `metadata()`.
#'
#' @export
#' @note Updated 2023-09-26.
#'
#' @return `ProteinToGene`.
setClass(
    Class = "ProteinToGene",
    contains = "DFrame"
)
setValidity(
    Class = "ProteinToGene",
    method = function(object) {
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
        ok <- validate(
            hasRows(object),
            all(complete.cases(object))
        )
        if (!isTRUE(ok)) {
            return(ok)
        }
        ok <- .validateMetadata(object)
        if (!isTRUE(ok)) {
            return(ok)
        }
        TRUE
    }
)



#' @inherit AcidGenerics::TxToGene description return title
#' @export
#' @note Updated 2023-09-26.
#'
#' @details
#' Contains a `DFrame` with `"txId"` and `"geneId"` columns.
#'
#' @section Genome metadata:
#'
#' We recommend slotting `"organism"`, `"genomeBuild"`, and `"release"` into
#' `metadata()`.
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
        ok <- validateClasses(
            object = object,
            expected = list(
                "txId" = "character",
                "geneId" = "character"
            )
        )
        if (!isTRUE(ok)) {
            return(ok)
        }
        ok <- validate(
            hasRows(object),
            all(complete.cases(object)),
            hasNoDuplicates(object[["txId"]])
        )
        if (!isTRUE(ok)) {
            return(ok)
        }
        ok <- .validateMetadata(object)
        if (!isTRUE(ok)) {
            return(ok)
        }
        TRUE
    }
)
