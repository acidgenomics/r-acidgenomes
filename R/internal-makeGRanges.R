#' Add broad class annotations
#'
#' @note Updated 2020-10-05.
#' @noRd
.addBroadClass <- function(object) {
    assert(is(object, "GRanges"))
    mcols(object)[["broadClass"]] <- Rle(.broadClass(object))
    object
}



#' Add Ensembl gene synonyms
#'
#' @note Updated 2020-10-05.
#' @noRd
.addGeneSynonyms <- function(object) {
    assert(is(object, "GRanges"))
    mcols <- mcols(object)
    assert(isSubset("geneID", colnames(mcols)))
    if (!any(grepl(pattern = "^ENS", x = mcols[["geneID"]]))) {
        return(object)
    }
    organism <- organism(object)
    if (!isSubset(organism, eval(formals(geneSynonyms)[["organism"]]))) {
        return(object)
    }
    alert(sprintf(
        "Adding gene synonyms to {.var %s} column.",
        "geneSynonyms"
    ))
    synonyms <- geneSynonyms(organism = organism, return = "DataFrame")
    assert(identical(c("geneID", "geneSynonyms"), colnames(synonyms)))
    mcols <- leftJoin(x = mcols, y = synonyms, by = "geneID")
    mcols(object) <- mcols
    object
}



#' Add the transcript version
#'
#' Append the transcript version to the identifier (e.g. ENST00000000233.10).
#'
#' @note Updated 2021-01-06.
#' @noRd
.addTxVersion <- function(object) {
    alert("Adding version to transcript identifiers.")
    mcolnames <- colnames(mcols(object))
    assert(
        is(object, "GRanges"),
        identical(metadata(object)[["level"]], "transcripts"),
        isSubset("transcriptID", mcolnames)
    )
    if (isSubset("transcriptIDVersion", mcolnames)) {
        ## `makeGRangesFromEnsembl()` output via ensembldb.
        id <- mcols(object)[["transcriptIDVersion"]]
    } else if (isSubset("transcriptVersion", mcolnames)) {
        ## `makeGRangesFromGFF()` output.
        id <- mcols(object)[["transcriptID"]]
        version <- mcols(object)[["transcriptVersion"]]
        id <- Rle(paste(id, version, sep = "."))
    } else {
        stop("Failed to locate transcript version metadata.")  # nolint
    }
    mcols(object)[["transcriptID"]] <- id
    ## Note that names are set by `.makeGRanges()`.
    object
}



#' Detect GRanges identifiers
#'
#' @note This intentionally prioritizes transcripts over genes.
#'
#' @note Updated 2020-10-05.
#' @noRd
.detectGRangesIDs <- function(object) {
    if (is(object, "GRangesList")) {
        object <- object[[1L]]
    }
    assert(is(object, "GRanges"))
    mcolnames <- colnames(mcols(object))
    if ("transcriptID" %in% mcolnames) {
        out <- "transcriptID"
    } else if ("transcript_id" %in% mcolnames) {
        out <- "transcript_id"
    } else if ("geneID" %in% mcolnames) {
        out <- "geneID"
    } else if ("gene_id" %in% mcolnames) {
        out <- "gene_id"
    } else {
        stop("Failed to detect ID column.")
    }
    out
}



## FIXME NEED TO IMPROVE CONSISTENCY Of METADATA RETURN.
## FIXME NEED TO DETECT ORGANISM AUTOMATICALLY.
## FIXME NEED TO SLOT ARGUMENTS INTO OBJECT.
## FIXME ORGANISM ALWAYS NEEDS TO BE DEFINED.
## FIXME DONT ALLOW GRANGESLIST RETURN HERE.

#' Make GRanges
#'
#' This is the main GRanges final return generator, used by
#' `makeGRangesFromEnsembl()` and `makeGRangesFromGFF()`.
#'
#' @note Updated 2020-10-05.
#' @noRd
.makeGRanges <- function(
    object,
    ignoreTxVersion = TRUE,
    broadClass = TRUE,
    synonyms = TRUE
) {
    assert(
        is(object, "GRanges"),
        hasNames(object),
        hasLength(object),
        isFlag(ignoreTxVersion),
        isFlag(broadClass),
        isFlag(synonyms)
    )
    length <- length(object)
    object <- .minimizeGRanges(object)
    assert(hasLength(object, n = length))
    object <- .standardizeGRanges(object)
    assert(hasLength(object, n = length))
    if (isFALSE(ignoreTxVersion)) {
        object <- .addTxVersion(object)
    }
    if (isTRUE(broadClass)) {
        object <- .addBroadClass(object)
    }
    if (isTRUE(synonyms)) {
        object <- .addGeneSynonyms(object)
    }
    ## Sort the metadata columns alphabetically.
    mcols(object) <-
        mcols(object)[, sort(colnames(mcols(object))), drop = FALSE]
    ## Prepare the metadata.
    ## Slot organism into metadata.
    object <- .slotOrganism(object)
    ## Ensure object contains prototype metadata.
    metadata(object) <- c(.prototypeMetadata, metadata(object))
    idCol <- .detectGRangesIDs(object)
    assert(isSubset(idCol, colnames(mcols(object))))
    names <- as.character(mcols(object)[[idCol]])
    assert(!any(is.na(names)))
    ## Inform the user if the object contains invalid names, showing offenders.
    invalid <- setdiff(names, make.names(names, unique = TRUE))
    if (hasLength(invalid)) {
        invalid <- sort(unique(invalid))
        alertWarning(sprintf(
            fmt = "%d invalid %s: {.var %s}.",
            length(invalid),
            ngettext(
                n = length(invalid),
                msg1 = "name",
                msg2 = "names"
            ),
            toInlineString(invalid, n = 10L)
        ))
    }
    rm(invalid)
    ## Split into GRangesList if object contains multiple ranges per feature.
    if (hasDuplicates(names)) {
        alertWarning(sprintf(
            fmt = paste(
                "{.var %s} contains multiple ranges per '%s'.",
                "Splitting into {.var %s}."
            ),
            "GRanges", idCol, "GRangesList"
        ))
        ## Metadata will get dropped during `split()` call; stash and reassign.
        metadata <- metadata(object)
        object <- split(x = object, f = as.factor(names))
        metadata(object) <- metadata
        rm(metadata)
    } else {
        names(object) <- names
    }
    ## Ensure the ranges are sorted by gene identifier.
    object <- object[sort(names(object))]
    ## Inform the user about the number of features returned.
    level <- match.arg(
        arg = metadata(object)[["level"]],
        choices = c("genes", "transcripts")
    )
    alertInfo(sprintf(
        "%d %s detected.",
        length(object),
        ngettext(
            n = length(object),
            msg1 = substr(level, 1L, nchar(level) - 1L),  # gene
            msg2 = level                                  # genes
        )
    ))

    ## FIXME ADD ASSERT FOR ORGANISM HERE.
    ## FIXME ADD ASSERT FOR CALL HERE TOO.

    object
}



#' Minimize GRanges
#'
#' This step drops extra columns in `mcols()` and applies run-length encoding,
#' to reduce memory overhead.
#'
#' Note that `removeNA()` call currently will error on complex columns. For
#' example, this will error on `CharacterList` columns returned from GENCODE
#' GFF3 file.
#'
#' @note Updated 2020-01-20.
#' @noRd
.minimizeGRanges <- function(object) {
    assert(is(object, "GRanges"))
    mcols <- mcols(object)
    ## Drop any complex S4 columns that aren't either atomic or list.
    keep <- bapply(
        X = mcols,
        FUN = function(x) {
            is.atomic(x) || is.list(x)
        }
    )
    mcols <- mcols[, keep, drop = FALSE]
    ## Ensure NA values are properly set, prior to `removeNA()` call.
    mcols <- sanitizeNA(mcols)
    ## Remove columns that are all `NA`. This step will remove all
    ## transcript-level columns from gene-level ranges.
    mcols <- removeNA(mcols)
    ## Apply run-length encoding on all atomic columns.
    mcols <- lapply(
        X = mcols,
        FUN = function(x) {
            if (isS4(x) || is(x, "AsIs") || !is.atomic(x)) {
                ## `I()` inhibits reinterpretation and returns `AsIs` class.
                ## This keeps complex columns (e.g. Entrez list) intact.
                ## Recommended in the `DataFrame` documentation.
                I(x)
            } else {
                ## Ensure factor levels get correctly reset, to save memory.
                if (is.factor(x)) {
                    x <- droplevels(x)
                }
                ## Use S4 run length encoding (Rle) for atomic metadata columns.
                ## Many of these elements are repetitive, and this makes
                ## operations faster.
                Rle(x)
            }
        }
    )
    ## `lapply()` returns as list, so we need to coerce back to DataFrame.
    mcols <- as(mcols, "DataFrame")
    mcols(object) <- mcols
    object
}



#' Standardize the GRanges into desired conventions
#'
#' Note that this step makes GRanges imported via `rtracklayer::import()`
#' incompatible with `GenomicFeatures::makeTxDbFromGRanges()`.
#'
#' @note Updated 2020-10-06.
#' @noRd
.standardizeGRanges <- function(object) {
    assert(is(object, "GRanges"))
    ## Standardize the metadata columns.
    mcols <- mcols(object)
    ## Use `transcript` prefix instead of `tx` consistently.
    colnames(mcols) <- gsub(
        pattern = "^tx_",
        replacement = "transcript_",
        x = colnames(mcols)
    )
    ## Ensure "ID" is always capitalized (e.g. "entrezid").
    colnames(mcols) <- gsub(
        pattern = "(.+)id$",
        replacement = "\\1ID",
        x = colnames(mcols)
    )
    ## Always return using camel case, even though GFF/GTF files use snake.
    colnames(mcols) <- camelCase(colnames(mcols))
    ## Always use `geneName` instead of `symbol`.
    ## Note that ensembldb output duplicates these.
    if (all(c("geneName", "symbol") %in% colnames(mcols))) {
        mcols[["symbol"]] <- NULL
    } else if ("symbol" %in% colnames(mcols)) {
        alert(
            "Renaming {.var symbol} to {.var geneName} in {.fun mcols}."
        )
        mcols[["geneName"]] <- mcols[["symbol"]]
        mcols[["symbol"]] <- NULL
    }
    ## Re-slot updated mcols back into object before calculating broad class
    ## biotype and/or assigning names.
    mcols(object) <- mcols
    ## Ensure the ranges are sorted by identifier.
    idCol <- .detectGRangesIDs(object)
    alert(sprintf("Arranging by {.var %s}.", idCol))
    names(object) <- mcols(object)[[idCol]]
    object <- object[sort(names(object))]
    object
}
