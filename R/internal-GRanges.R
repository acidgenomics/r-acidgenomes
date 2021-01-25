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
#' @note Updated 2021-01-18.
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
    } else if (x[["biotype"]] == "protein_coding") {
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



#' Add broad class annotations
#'
#' @note Updated 2021-01-25.
#' @noRd
.addBroadClass <- function(object) {
    assert(
        is(object, "GRanges"),
        hasNoDuplicates(names(object)),
        identical(
            x = colnames(mcols(object)),
            y = camelCase(colnames(mcols(object)), strict = TRUE)
        ),
        !isSubset("broadClass", colnames(mcols(object)))
    )
    df <- as.data.frame(object)
    colnames(df) <- gsub(
        pattern = "^transcript",
        replacement = "tx",
        x = colnames(df)
    )
    ## Biotypes ----------------------------------------------------------------
    ## Prioritizing transcript biotype over gene, if defined. This only applies
    ## for transcript-level GRanges. For gene-level GRanges, the gene biotypes
    ## will be used, as expected.
    if ("txBiotype" %in% colnames(df)) {
        biotypeCol <- "txBiotype"
        biotypeData <- df[[biotypeCol]]
    } else if ("geneBiotype" %in% colnames(df)) {
        biotypeCol <- "geneBiotype"
        biotypeData <- df[[biotypeCol]]
    } else {
        return(NA_character_)
    }
    ## Gene names --------------------------------------------------------------
    if ("geneName" %in% colnames(df)) {
        geneNameCol <- "geneName"
        geneNameData <- df[[geneNameCol]]
    } else {
        geneNameCol <- NULL
        geneNameData <- NA_character_
    }
    ## Seqnames ----------------------------------------------------------------
    ## This refers to the chromosome name.
    ## Note that data frame coercion will define `seqnames` column from the
    ## `GRanges` object (see above).
    if ("seqnames" %in% colnames(df)) {
        seqnamesCol <- "seqnames"
        seqnamesData <- df[[seqnamesCol]]
    } else {
        seqnamesCol <- NULL
        seqnamesData <- NA_character_
    }
    ## Apply broad class -------------------------------------------------------
    ## Note that this method doesn't seem to work right with DataFrame class.
    df <- data.frame(
        "biotype" = biotypeData,
        "chromosome" = seqnamesData,
        "geneName" = geneNameData,
        stringsAsFactors = TRUE
    )
    ## Consider adding BiocParallel support here for improved speed.
    mcols(object)[["broadClass"]] <-
        Rle(as.factor(apply(X = df, MARGIN = 1L, FUN = .applyBroadClass)))
    validObject(object)
    object
}



#' Add Ensembl gene synonyms
#'
#' @note Updated 2021-01-18.
#' @noRd
.addGeneSynonyms <- function(object) {
    assert(is(object, "GRanges"))
    mcols <- mcols(object)
    assert(
        isSubset("geneId", colnames(mcols)),
        allAreMatchingRegex(x = mcols[["geneId"]], pattern = "^ENS"),
    )
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



#' Standardize the GRanges mcols into desired naming conventions
#'
#' @details
#' Always return using camel case, even though GFF/GTF files use snake.
#'
#' Note that this step makes GRanges imported via `rtracklayer::import()`
#' incompatible with `GenomicFeatures::makeTxDbFromGRanges()` parser, so be
#' sure to call that function prior to attempting to run this step.
#'
#' @note Updated 2021-01-25.
#' @noRd
.standardizeGRanges <- function(object) {
    assert(is(object, "GRanges"))
    mcols <- mcols(object)
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
    ## Always prefer use "geneName" instead of "symbol".
    ## Note that ensembldb output duplicates these by default.
    if (isSubset(c("geneName", "symbol"), colnames(mcols))) {
        mcols[["symbol"]] <- NULL
    } else if (isSubset("symbol", colnames(mcols))) {
        mcols[["geneName"]] <- mcols[["symbol"]]
        mcols[["symbol"]] <- NULL
    }
    ## Remove any uninformative blacklisted columns.
    blacklistCols <- c(
        ## e.g. Ensembl GFF. Use "gene_biotype", "tx_biotype" instead.
        "biotype",
        ## e.g. Ensembl GFF: "havana_homo_sapiens". Not informative.
        "logic_name",
        "type"
        ## FIXME Other values to consider:
        ## "end_range",
        ## "exception",
        ## "gbkey",
        ## "partial",
        ## "pseudo",
        ## "start_range",
        ## "transl_except"

    )
    keep <- !colnames(mcols) %in% blacklistCols
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
#' @note Updated 2021-01-25.
#' @noRd
.makeGRanges <- function(
    object,
    ignoreVersion = TRUE,
    synonyms = FALSE,
    ## Internal-only arguments:
    broadClass = TRUE
) {
    assert(
        is(object, "GRanges"),
        hasLength(object),
        isFlag(ignoreVersion),
        isFlag(synonyms),
        isFlag(broadClass),
        isString(metadata(object)[["level"]]),
        isString(metadata(object)[["provider"]])
    )
    level <- match.arg(
        arg = metadata(object)[["level"]],
        choices = .grangesLevels
    )
    object <- .minimizeGRanges(object)
    object <- .standardizeGRanges(object)
    if (isFALSE(ignoreVersion)) {
        object <- .addGeneVersion(object)
        object <- .addTxVersion(object)
    }
    if (isTRUE(broadClass)) {
        object <- .addBroadClass(object)
    }
    if (isTRUE(synonyms)) {
        object <- .addGeneSynonyms(object)
    }
    ## Define the names by desired identifier column.
    idCol <- .matchGRangesNamesColumn(object)
    assert(isSubset(idCol, colnames(mcols(object))))
    alert(sprintf("Defining names by {.var %s} column.", idCol))
    names <- as.character(mcols(object)[[idCol]])
    assert(hasNoDuplicates(names), !any(is.na(names)))
    ## Inform the user if the object contains invalid names, showing offenders.
    ## This can happen with RefSeq genes, WormBase transcripts, but should be
    ## clean for Ensembl and GENCODE.
    invalidNames <- setdiff(names, make.names(names, unique = TRUE))
    if (hasLength(invalidNames)) {
        invalidNames <- sort(unique(invalidNames))
        alertWarning(sprintf(
            fmt = "%d invalid %s: %s.",
            length(invalidNames),
            ngettext(
                n = length(invalidNames),
                msg1 = "name",
                msg2 = "names"
            ),
            toInlineString(invalidNames, n = 10L)
        ))
    }
    names(object) <- names
    ## Ensure the ranges are sorted by genomic position.
    object <- sort(object)
    assert(isFALSE(is.unsorted(object)))
    ## Inform the user about the number of features returned.
    alertInfo(sprintf(
        "%d %s detected.",
        length(object),
        ngettext(
            n = length(object),
            msg1 = substr(level, 1L, nchar(level) - 1L),  # gene
            msg2 = level                                  # genes
        )
    ))
    ## Sort the metadata columns alphabetically.
    mcols(object) <-
        mcols(object)[, sort(colnames(mcols(object))), drop = FALSE]
    ## Ensure metadata elements are all sorted alphabetically.
    metadata(object)[["ignoreVersion"]] <- ignoreVersion
    metadata(object) <- metadata(object)[sort(names(metadata(object)))]
    ## Run final assert checks before returning.
    validObject(object)
    provider <- metadata(object)[["provider"]]
    assert(isString(provider))
    new(Class = upperCamelCase(paste(provider, level)), object)
}
