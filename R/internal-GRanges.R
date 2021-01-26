## Metadata modification =======================================================
#' Apply broad class definitions
#'
#' This function is intended to work rowwise on the GRanges mcols.
#'
#' @section Mitochondrial genes:
#' Mitochondrial gene matching depends on the genome.
#' This is important in particular for single-cell RNA-seq.
#'
#' - *H. sapiens*: "MT-" gene name.
#' - *M. musculus*: "mt-" gene name.
#' - *D. melanogaster*: "mt:" gene name.
#' - *C. elegans*: Can't match by gene name. Match by "MtDNA" chromosome.
#'   Alternatively, can match using "MTCE" sequence name (parent clone).
#'   https://www.wormbase.org/species/c_elegans/clone/MTCE
#'
#' Note that this might not be perfect for other genomes, so consider
#' atttempting to improve support here in a future update.
#'
#' @seealso Can use `dplyr::case_when()` instead, which allows for a rowwise
#'   vectorized if/else call stack.
#'
#' @note Can return `NA_character_` here instead. Keeping this as "other", to
#'   main consistency with previous data sets. Also note that `NA` can behave
#'   inconsistently in plotting engines.
#' @note Updated 2021-01-25.
#'
#' @author Rory Kirchner, Michael Steinbaugh
#' @noRd
#'
#' @param x `list`.
#'   List returned via apply call using `MARGIN = 1`.
#'
#' @return `character(1)`.
.applyBroadClass <- function(x) {
    if (
        isTRUE(grepl(pattern = "^MT",
            x = x[["chromosome"]],
            ignore.case = TRUE
        )) ||
        isTRUE(grepl(
            pattern = "^mt[\\:\\-]",
            x = x[["geneName"]],
            ignore.case = TRUE
        ))
    ) {
        "mito"
    } else if (identical(
        x = x[["biotype"]],
        y = "protein_coding"
    )) {
        "coding"
    } else if (
        x[["biotype"]] %in% c(
            "known_ncrna",
            "lincRNA",
            "non_coding"
        )
    ) {
        "noncoding"
    } else if (
        isTRUE(grepl(
            pattern = "pseudo",
            x = x[["biotype"]],
            ignore.case = TRUE
        ))
    ) {
        "pseudo"
    } else if (
        x[["biotype"]] %in% c(
            "miRNA",
            "misc_RNA",
            "ribozyme",
            "rRNA",
            "scaRNA",
            "scRNA",
            "snoRNA",
            "snRNA",
            "sRNA"
        )
    ) {
        "small"
    } else if (
        x[["biotype"]] %in% c(
            "non_stop_decay",
            "nonsense_mediated_decay"
        )
    ) {
        "decaying"
    } else if (
        isTRUE(grepl(
            pattern = "^ig_",
            x = x[["biotype"]],
            ignore.case = TRUE
        ))
    ) {
        "ig"  ## immunoglobulin
    } else if (
        isTRUE(grepl(
            pattern = "^tr_",
            x = x[["biotype"]],
            ignore.case = TRUE
        ))
    ) {
        "tcr"  ## T cell receptor
    } else {
        "other"
    }
}



## FIXME Make this tighter and ONLY return if biotype column is defined.

#' Add broad class annotations
#'
#' @note Updated 2021-01-25.
#' @noRd
.addBroadClass <- function(object) {
    assert(
        is(object, "GRanges"),
        identical(
            x = colnames(mcols(object)),
            y = camelCase(colnames(mcols(object)), strict = TRUE)
        ),
        !isSubset("broadClass", colnames(mcols(object))),
        allAreNotMatchingRegex(
            pattern = "^transcript",
            x = colnames(mcols(object))
        )
    )
    df <- as.data.frame(object)
    ## Biotypes. Prioritizing transcript biotype over gene, if defined. This
    ## only applies for transcript-level GRanges. For gene-level GRanges, the
    ## gene biotypes will be used, as expected.
    if ("txBiotype" %in% colnames(df)) {
        biotypeCol <- "txBiotype"
        biotypeData <- df[[biotypeCol]]
    } else if ("geneBiotype" %in% colnames(df)) {
        biotypeCol <- "geneBiotype"
        biotypeData <- df[[biotypeCol]]
    } else {
        biotypeCol <- NULL
        biotypeData <- NA_character_
    }
    ## Gene names.
    if ("geneName" %in% colnames(df)) {
        geneNameCol <- "geneName"
        geneNameData <- df[[geneNameCol]]
    } else {
        geneNameCol <- NULL
        geneNameData <- NA_character_
    }
    ## Seqnames. This refers to the chromosome name. Note that data frame
    ## coercion will define `seqnames` column from the `GRanges` object
    ## (see above).
    if ("seqnames" %in% colnames(df)) {
        seqnamesCol <- "seqnames"
        seqnamesData <- df[[seqnamesCol]]
    } else {
        seqnamesCol <- NULL
        seqnamesData <- NA_character_
    }
    ## Apply broad class. Note that this method doesn't seem to work right with
    ## DataFrame class.
    df <- data.frame(
        "biotype" = biotypeData,
        "chromosome" = seqnamesData,
        "geneName" = geneNameData,
        stringsAsFactors = TRUE
    )
    ## NOTE Consider using BiocParallel here for improved speed.
    x <- apply(X = df, MARGIN = 1L, FUN = .applyBroadClass)
    x <- as.factor(x)
    if (all(x == "other")) {
        return(object)
    }
    mcols(object)[["broadClass"]] <- Rle(x)
    object
}



## FIXME ONLY ALLOW THIS FOR ENSEMBL AND GENCODE.

#' Add Ensembl gene synonyms
#'
#' @note Updated 2021-01-26.
#' @noRd
#'
#' @details Currently supported only for Ensembl and GENCODE genomes.
.addGeneSynonyms <- function(object) {
    assert(
        is(object, "GRanges"),
        isString(metadata(object)[["provider"]]),
        isSubset(
            x = metadata(object)[["provider"]],
            y = c("Ensembl", "GENCODE")
        ),
        isSubset("geneId", colnames(mcols(object))),
        allAreMatchingRegex(x = mcols(object)[["geneId"]], pattern = "^ENS")
    )
    mcols <- mcols(object)
    organism <- match.arg(
        arg = organism(object),
        choices = eval(formals(geneSynonyms)[["organism"]])
    )
    alert(sprintf(
        "Adding gene synonyms for {.var %s} to {.var %s} column.",
        organism, "geneSynonyms"
    ))
    synonyms <- geneSynonyms(organism = organism, return = "DataFrame")
    assert(identical(c("geneId", "geneSynonyms"), colnames(synonyms)))
    mcols <- leftJoin(x = mcols, y = synonyms, by = "geneId")
    mcols(object) <- mcols
    object
}



## FIXME Should this error out for unsupported genome?

#' Add the gene identifier version
#'
#' Append the gene version to the identifier (e.g. ENSG00000000003.15).
#'
#' @note Updated 2021-01-21.
#' @noRd
.addGeneVersion <- function(object) {
    assert(is(object, "GRanges"))
    if (!isSubset("geneId", colnames(mcols(object)))) {
        return(object)
    }
    alert("Including version in gene identifiers.")
    if (isSubset("geneIdVersion", colnames(mcols(object)))) {
        id <- mcols(object)[["geneIdVersion"]]
    } else if (isSubset("geneVersion", colnames(mcols(object)))) {
        id <- Rle(paste(
            mcols(object)[["geneId"]],
            mcols(object)[["geneVersion"]],
            sep = "."
        ))
    } else {
        stop("Failed to locate gene identifier version.")  # nocov
    }
    mcols(object)[["geneIdNoVersion"]] <- mcols(object)[["geneId"]]
    mcols(object)[["geneId"]] <- id
    object
}



## FIXME Should this error out for unsupported genome?

#' Add the transcript identifier version
#'
#' Append the transcript version to the identifier (e.g. ENST00000000233.10).
#'
#' @note Updated 2021-01-21.
#' @noRd
.addTxVersion <- function(object) {
    assert(is(object, "GRanges"))
    if (!isSubset("txId", colnames(mcols(object)))) {
        return(object)
    }
    alert("Including version in transcript identifiers.")
    if (isSubset("txIdVersion", colnames(mcols(object)))) {
        id <- mcols(object)[["txIdVersion"]]
    } else if (isSubset("txVersion", colnames(mcols(object)))) {
        id <- Rle(paste(
            mcols(object)[["txId"]],
            mcols(object)[["txVersion"]],
            sep = "."
        ))
    } else {
        stop("Failed to locate transcript identifier version.")  # nocov
    }
    mcols(object)[["txIdNoVersion"]] <- mcols(object)[["txId"]]
    mcols(object)[["txId"]] <- id
    object
}



## Standardization =============================================================
#' Match the identifier column in GRanges to use for names.
#'
#' @note Updated 2021-01-20.
#' @noRd
.matchGRangesNamesColumn <- function(object) {
    assert(
        is(object, "GRanges"),
        isString(metadata(object)[["level"]])
    )
    level <- match.arg(
        arg = metadata(object)[["level"]],
        choices = c(
            "cds",
            "exons",
            "genes",
            "transcripts"
        )
    )
    x <- camelCase(colnames(mcols(object)), strict = TRUE)
    table <- switch(
        EXPR = level,
        "cds" = "cdsId",
        "exons" = "exonId",
        "genes" = "geneId",
        "transcripts" = "txId"
    )
    idx <- which(x %in% table)
    x <- colnames(mcols(object))[idx]
    assert(isString(x))
    x
}



#' Minimize GRanges mcols
#'
#' This step sanitizes NA values, applies run-length encoding (to reduce memory
#' overhead), and trims any invalid ranges.
#'
#' @note Updated 2021-01-25.
#' @noRd
.minimizeGRanges <- function(object) {
    assert(is(object, "GRanges"))
    ## This trimming step was added to handle GRanges from Ensembl 102, which
    ## won't return valid otherwise from ensembldb.
    length <- length(object)
    object <- trim(object)
    assert(hasLength(object, n = length))
    mcols <- mcols(object)
    mcolsList <- lapply(
        X = mcols,
        FUN = function(x) {
            if (isS4(x) || is(x, "AsIs") || !is.atomic(x)) {
                I(x)
            } else {
                x <- sanitizeNA(x)
                if (all(is.na(x))) {
                    return(NULL)
                }
                if (is.factor(x)) {
                    x <- droplevels(x)
                }
                x <- Rle(x)
                x
            }
        }
    )
    mcolsList <- Filter(f = Negate(is.null), x = mcolsList)
    mcols <- as(mcolsList, "DataFrame")
    mcols(object) <- mcols
    object
}



## FIXME NEED TO AUTOMATICALLY NUKE CAPITZLIED COLUMNS.

## FIXME NEED TO ADD WHITELIST COLUMNS.
##       Dbxref, db_xref from RefSeq...need to standardize the convention here.
##       Then safe to nuke all capital columns.

#' Standardize the GRanges mcols into desired naming conventions
#'
#' @details
#' Always return using camel case, even though GFF/GTF files use snake.
#'
#' Note that this step makes GRanges imported via `rtracklayer::import()`
#' incompatible with `GenomicFeatures::makeTxDbFromGRanges()` parser, so be
#' sure to call that function prior to attempting to run this step.
#'
#' @note Updated 2021-01-26.
#' @noRd
.standardizeGRanges <- function(object) {
    assert(is(object, "GRanges"))
    mcols <- mcols(object)


    ## FIXME AUTOMATICALLY REMOVE CAPITALIZED COLUMNS, EXCEPT FOR A WHITELIST.



    ## Changed to strict format here in v0.2.0 release.
    ## This results in returning "Id" identifier suffix instead of "ID".
    colnames(mcols) <- camelCase(colnames(mcols), strict = TRUE)
    ## Ensure "tx" prefix is used consistently instead of "transcript".
    ## This convention was changed in v0.2.0 release.
    colnames(mcols) <- gsub(
        pattern = "^transcript",
        replacement = "tx",
        x = colnames(mcols),
        ignore.case = FALSE
    )
    ## Ensure "Id" is always capitalized (e.g. "entrezid" to "entrezId").
    colnames(mcols) <-
        gsub(
            pattern = "(.+)id$",
            replacement = "\\1Id",
            x = colnames(mcols),
            ignore.case = FALSE
    )
    ## Remove any uninformative blacklisted columns.
    blacklistCols <- c(
        ## e.g. Ensembl GFF. Use "(gene|tx)_biotype" instead.
        "biotype",
        ## e.g. RefSeq GFF.
        "gbkey",
        ## e.g. Ensembl GFF: "havana_homo_sapiens". Not informative.
        "logic_name",
        "type"
        ## FIXME Other values to consider:
        ## "end_range",
        ## "exception",
        ## "partial",
        ## "pseudo",
        ## "start_range",
        ## "transl_except"
    )
    keep <- !colnames(mcols) %in% blacklistCols
    ## Always prefer use of "geneName" instead of "geneSymbol" or "symbol".
    ## Note that ensembldb output "symbol" duplicate by default.
    if (isSubset(c("geneName", "symbol"), colnames(mcols))) {
        mcols[["symbol"]] <- NULL
    } else if (
        isSubset("symbol", colnames(mcols)) &&
        !isSubset("geneName", colnames(mcols))
    ) {
        colnames(mcols)[colnames(mcols) == "symbol"] <- "geneName"
    } else if (
        isSubset("geneSymbol", colnames(mcols)) &&
        !isSubset("geneName", colnames(mcols))
    ) {
        ## e.g. FlyBase GTF.
        colnames(mcols)[colnames(mcols) == "geneSymbol"] <- "geneName"
    }
    ## Always prefer use of "txName" instead of "txSymbol".
    if (
        isSubset("txSymbol", colnames(mcols)) &&
        !isSubset("txName", colnames(mcols))
    ) {
        ## e.g. FlyBase GTF.
        colnames(mcols)[colnames(mcols) == "txSymbol"] <- "txName"
    }
    ## Add geneName column if missing.
    if (
        isSubset("geneId", colnames(mcols)) &&
        !isSubset("geneName", colnames(mcols))
    ) {
        mcols[["geneName"]] <- mcols[["geneId"]]
    }
    ## Add txName column if missing.
    if (
        isSubset("txId", colnames(mcols)) &&
        !isSubset("txName", colnames(mcols))
    ) {
        mcols[["txName"]] <- mcols[["txId"]]
    }
    ## Always prefer use of "geneBiotype" instead of "geneType".
    if (
        isSubset("geneType", colnames(mcols)) &&
        !isSubset("geneBiotype", colnames(mcols))
    ) {
        ## e.g. GENCODE GFF.
        colnames(mcols)[colnames(mcols) == "geneType"] <- "geneBiotype"
    }
    if (
        isSubset("txType", colnames(mcols)) &&
        !isSubset("txBiotype", colnames(mcols))
    ) {
        ## e.g. GENCODE GFF.
        colnames(mcols)[colnames(mcols) == "txType"] <- "txBiotype"
    }
    mcols <- mcols[keep]
    mcols(object) <- mcols
    object
}



## Main generator ==============================================================
#' Make GRanges
#'
#' This is the main GRanges final return generator, used by
#' `makeGRangesFromEnsembl()` and `makeGRangesFromGFF()`.
#'
#' @note Updated 2021-01-26.
#' @noRd
.makeGRanges <- function(
    object,
    ignoreVersion = TRUE,
    synonyms = FALSE
) {
    assert(
        is(object, "GRanges"),
        hasLength(object),
        isFlag(ignoreVersion),
        isFlag(synonyms),
        isString(metadata(object)[["level"]]),
        isString(metadata(object)[["provider"]])
    )
    level <- match.arg(
        arg = metadata(object)[["level"]],
        choices = .grangesLevels
    )
    provider <- metadata(object)[["provider"]]
    object <- .minimizeGRanges(object)
    object <- .standardizeGRanges(object)
    if (isFALSE(ignoreVersion)) {
        object <- .addGeneVersion(object)
        object <- .addTxVersion(object)
    }
    object <- .addBroadClass(object)
    if (isTRUE(synonyms)) {
        object <- .addGeneSynonyms(object)
    }
    idCol <- .matchGRangesNamesColumn(object)
    assert(isSubset(idCol, colnames(mcols(object))))
    alert(sprintf("Defining names by {.var %s} column.", idCol))
    if (hasDuplicates(mcols(object)[[idCol]])) {
        alertInfo(sprintf(
            fmt = paste(
                "{.var %s} contains multiple ranges per {.var %s}.",
                "Splitting into {.var %s}."
            ),
            "GRanges", idCol, "GRangesList"
        ))
        ## Metadata will get dropped during `split()` call; stash and reassign.
        meta <- metadata(object)
        object <- split(x = object, f = as.factor(mcols(object)[[idCol]]))
        metadata(object) <- meta
    } else {
        names <- as.character(mcols(object)[[idCol]])
        assert(hasNoDuplicates(names), !any(is.na(names)))
        names(object) <- names
        object <- sort(object)
        assert(isFALSE(is.unsorted(object)))
    }
    ## Inform the user about the number of features returned.
    ## > alertInfo(sprintf(
    ## >     "%d %s detected.",
    ## >     length(object),
    ## >     ngettext(
    ## >         n = length(object),
    ## >         msg1 = substr(level, 1L, nchar(level) - 1L),  # gene
    ## >         msg2 = level                                  # genes
    ## >     )
    ## > ))
    ## Sort the mcols alphabetically.
    mcols(object) <-
        mcols(object)[, sort(colnames(mcols(object))), drop = FALSE]
    ## Ensure metadata elements are all sorted alphabetically.
    metadata(object)[["ignoreVersion"]] <- ignoreVersion
    metadata(object) <- metadata(object)[sort(names(metadata(object)))]
    ## Run final assert checks before returning.
    validObject(object)
    class <- upperCamelCase(
        object = paste(provider, level),
        strict = FALSE
    )
    if (isClass(Class = class)) {
        out <- new(Class = class, object)
    } else {
        ## This is used to return CDS and exons from TxDb.
        out <- object
    }
    out
}
