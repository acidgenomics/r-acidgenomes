## Metadata modification =======================================================

## nolint start

#' Apply broad class definitions
#'
#' This function is intended to work rowwise on the GRanges mcols.
#'
#' @section Mitochondrial genes:
#'
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
#' @section Ensembl biotypes:
#'
#' See [biotypes guide](https://m.ensembl.org/info/genome/genebuild/biotypes.html)
#' for details.
#'
#' LRG: Locus Reference Genomic sequence.
#' Refer to the [LRG website](https://www.lrg-sequence.org/) for details.
#'
#' TEC (To be Experimentally Confirmed): Regions with EST clusters that have
#' polyA features that could indicate the presence of protein coding genes.
#' These require experimental validation, either by 5' RACE or RT-PCR to extend
#' the transcripts, or by confirming expression of the putatively-encoded
#' peptide with specific antibodies.
#'
#' @seealso Can use `dplyr::case_when()` instead, which allows for a rowwise
#'   vectorized if/else call stack.
#'
#' @note Can return `NA_character_` here instead. Keeping this as "other", to
#'   main consistency with previous data sets. Also note that `NA` can behave
#'   inconsistently in plotting engines.
#' @note Updated 2021-01-30.
#'
#' @author Rory Kirchner, Michael Steinbaugh
#' @noRd
#'
#' @param x `list`.
#'   List returned via apply call using `MARGIN = 1`.
#'
#' @return `character(1)`.

## nolint end

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
            "lncRNA",
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
#' @note Updated 2021-01-30.
#' @noRd
.addBroadClass <- function(object) {
    assert(
        is(object, "GRanges"),
        identical(
            x = names(mcols(object)),
            y = camelCase(names(mcols(object)), strict = TRUE)
        ),
        !isSubset("broadClass", names(mcols(object))),
        allAreNotMatchingRegex(
            pattern = "^transcript",
            x = names(mcols(object))
        )
    )
    df <- as.data.frame(object)
    ## Biotypes. Prioritizing gene over transcript biotype, if defined. This
    ## only applies for transcript-level GRanges.
    if ("geneBiotype" %in% names(df)) {
        biotypeCol <- "geneBiotype"
        biotypeData <- df[[biotypeCol]]
    } else if ("txBiotype" %in% names(df)) {
        biotypeCol <- "txBiotype"
        biotypeData <- df[[biotypeCol]]
    } else {
        biotypeCol <- NULL
        biotypeData <- NA_character_
    }
    ## Gene names.
    if ("geneName" %in% names(df)) {
        geneNameCol <- "geneName"
        geneNameData <- df[[geneNameCol]]
    } else {
        geneNameCol <- NULL
        geneNameData <- NA_character_
    }
    ## Seqnames. This refers to the chromosome name. Note that data frame
    ## coercion will define `seqnames` column from the `GRanges` object
    ## (see above).
    if ("seqnames" %in% names(df)) {
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
        row.names = names(object),
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



#' Add Ensembl gene synonyms
#'
#' @note Updated 2021-01-27.
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
        isSubset("geneId", names(mcols(object))),
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
    assert(identical(c("geneId", "geneSynonyms"), names(synonyms)))
    mcols <- leftJoin(x = mcols, y = synonyms, by = "geneId")
    mcols(object) <- mcols
    object
}



## Identifier versions =========================================================
#' Include identifier version in primary identifier
#'
#' @note Updated 2021-01-30.
#' @noRd
.includeVersion <-
    function(
        object,
        idCol,
        idVersionCol,
        idNoVersionCol,
        quiet = TRUE
    ) {
        assert(
            is(object, "GRanges"),
            isString(idCol),
            isString(idVersionCol),
            isString(idNoVersionCol),
            isFlag(quiet)
        )
        ## Early return for genomes without identifier versions (e.g. FlyBase).
        if (!isSubset(
            x = c(idCol, idVersionCol),
            y = names(mcols(object))
        )) {
            return(object)
        }
        assert(areDisjointSets(idNoVersionCol, names(mcols(object))))
        if (isFALSE(quiet)) {
            alert(sprintf(
                "Including version in {.var %s} from {.var %s}.",
                idCol, idVersionCol
            ))
            alertInfo(sprintf(
                "Unversioned identifiers are in {.var %s}.",
                idNoVersionCol
            ))
        }
        id <- mcols(object)[[idVersionCol]]
        mcols(object)[[idNoVersionCol]] <- mcols(object)[[idCol]]
        mcols(object)[[idCol]] <- id
        object
    }



#' Include the gene identifier version
#'
#' Append the gene version to the identifier (e.g. ENSG00000000003.15).
#'
#' @note Updated 2021-01-27.
#' @noRd
.includeGeneVersion <- function(object) {
    .includeVersion(
        object = object,
        idCol = "geneId",
        idVersionCol = "geneIdVersion",
        idNoVersionCol = "geneIdNoVersion"
    )
}



#' Include the transcript identifier version
#'
#' Append the transcript version to the identifier (e.g. ENST00000000233.10).
#'
#' @note Updated 2021-01-27.
#' @noRd
.includeTxVersion <- function(object) {
    .includeVersion(
        object = object,
        idCol = "txId",
        idVersionCol = "txIdVersion",
        idNoVersionCol = "txIdNoVersion"
    )
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
    x <- camelCase(names(mcols(object)), strict = TRUE)
    table <- switch(
        EXPR = level,
        "cds" = "cdsId",
        "exons" = "exonId",
        "genes" = "geneId",
        "transcripts" = "txId"
    )
    idx <- which(x %in% table)
    x <- names(mcols(object))[idx]
    assert(isString(x))
    x
}



#' Minimize GRanges mcols
#'
#' This step sanitizes NA values, applies run-length encoding (to reduce memory
#' overhead), and trims any invalid ranges.
#'
#' @note Updated 2021-01-26.
#' @noRd
.minimizeMcols <- function(object) {
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
#' @note Updated 2021-01-26.
#' @noRd
.standardizeMcols <- function(object) {
    assert(is(object, "GRanges"))
    mcols <- mcols(object)
    ## Remove any columns beginning with a capital letter, which are used in
    ## GFF3 files.
    keep <- !grepl(pattern = "^[A-Z]", x = names(mcols))
    mcols <- mcols[keep]
    ## Changed to strict format here in v0.2.0 release. This results in
    ## returning "Id" identifier suffix instead of "ID".
    names(mcols) <- camelCase(names(mcols), strict = TRUE)
    ## Ensure "tx" prefix is used consistently instead of "transcript".
    ## This convention was changed in v0.2.0 release.
    names(mcols) <- gsub(
        pattern = "^transcript",
        replacement = "tx",
        x = names(mcols),
        ignore.case = FALSE
    )
    ## Ensure "Id" is always capitalized (e.g. "entrezid" to "entrezId").
    names(mcols) <-
        gsub(
            pattern = "(.+)id$",
            replacement = "\\1Id",
            x = names(mcols),
            ignore.case = FALSE
    )
    ## Always prefer use of "geneName" instead of "geneSymbol" or "symbol".
    ## Note that ensembldb output "symbol" duplicate by default.
    if (isSubset(c("geneName", "symbol"), names(mcols))) {
        mcols[["symbol"]] <- NULL
    } else if (
        isSubset("symbol", names(mcols)) &&
        !isSubset("geneName", names(mcols))
    ) {
        names(mcols)[names(mcols) == "symbol"] <- "geneName"
    } else if (
        isSubset("geneSymbol", names(mcols)) &&
        !isSubset("geneName", names(mcols))
    ) {
        ## e.g. FlyBase GTF.
        names(mcols)[names(mcols) == "geneSymbol"] <- "geneName"
    }
    ## Always prefer use of "txName" instead of "txSymbol".
    if (
        isSubset("txSymbol", names(mcols)) &&
        !isSubset("txName", names(mcols))
    ) {
        ## e.g. FlyBase GTF.
        names(mcols)[names(mcols) == "txSymbol"] <- "txName"
    }
    ## Add geneName column if missing.
    if (
        isSubset("geneId", names(mcols)) &&
        !isSubset("geneName", names(mcols))
    ) {
        mcols[["geneName"]] <- mcols[["geneId"]]
    }
    ## Add txName column if missing.
    if (
        isSubset("txId", names(mcols)) &&
        !isSubset("txName", names(mcols))
    ) {
        mcols[["txName"]] <- mcols[["txId"]]
    }
    ## Always prefer use of "geneBiotype" instead of "geneType" or "biotype".
    if (
        isSubset("geneType", names(mcols)) &&
        !isSubset("geneBiotype", names(mcols))
    ) {
        ## e.g. GENCODE GFF.
        names(mcols)[names(mcols) == "geneType"] <- "geneBiotype"
    } else if (
        isSubset("biotype", names(mcols)) &&
        !isSubset("geneBiotype", names(mcols))
    ) {
        names(mcols)[names(mcols) == "biotype"] <- "geneBiotype"
    }
    if (
        isSubset("txType", names(mcols)) &&
        !isSubset("txBiotype", names(mcols))
    ) {
        ## e.g. GENCODE GFF.
        names(mcols)[names(mcols) == "txType"] <- "txBiotype"
    }
    mcols(object) <- mcols
    object
}



## Main generator ==============================================================
#' Make GRanges
#'
#' This is the main GRanges final return generator, used by
#' `makeGRangesFromEnsembl()` and `makeGRangesFromGFF()`.
#'
#' @note Updated 2021-01-29.
#' @noRd
.makeGRanges <- function(
    object,
    ignoreVersion = FALSE,
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
    ## Consider not allowing user to ignore identifier versions for unsupported
    ## providers. Not failing here at the moment, to provide legacy support for
    ## bcbio R packages.
    ## > if (
    ## >     isTRUE(ignoreVersion) &&
    ## >     !isSubset(provider, c("Ensembl", "GENCODE"))
    ## > ) {
    ## >     stop(sprintf(
    ## >         paste(
    ## >             "Identifier version modification with '%s' flag",
    ## >             "is not supported for %s genomes."
    ## >         ),
    ## >         "ignoreVersion = TRUE",
    ## >         provider
    ## >     ))
    ## > }
    object <- .minimizeMcols(object)
    object <- .standardizeMcols(object)
    if (isFALSE(ignoreVersion)) {
        object <- .includeGeneVersion(object)
        object <- .includeTxVersion(object)
    }
    object <- .addBroadClass(object)
    if (isTRUE(synonyms)) {
        object <- .addGeneSynonyms(object)
    }
    idCol <- .matchGRangesNamesColumn(object)
    assert(isSubset(idCol, names(mcols(object))))
    alert(sprintf("Defining names by {.var %s} column.", idCol))
    ## Sort the ranges by genomic location.
    ## Previously we sorted by the identifier column, until v0.2.0.
    object <- sort(object)
    ## Sort the mcols alphabetically.
    mcols(object) <-
        mcols(object)[, sort(names(mcols(object))), drop = FALSE]
    ## Ensure metadata elements are all sorted alphabetically.
    metadata(object) <- append(
        x = metadata(object),
        values = list(
            "acidGenomes" = .version,
            "date" = Sys.Date(),
            "ignoreVersion" = ignoreVersion,
            "synonyms" = synonyms
        )
    )
    metadata(object) <- metadata(object)[sort(names(metadata(object)))]
    if (
        isSubset(provider, "RefSeq") ||
        hasDuplicates(mcols(object)[[idCol]])
    ) {
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
        assert(
            hasNoDuplicates(names),
            !any(is.na(names)),
            ## This check fails for split GRangesList (see above).
            isFALSE(is.unsorted(object))
        )
        names(object) <- names
    }
    ## Run final assert checks before returning.
    validObject(object)
    if (isSubset(level, c("genes", "transcripts"))) {
        class <- upperCamelCase(
            object = paste(
                switch(
                    EXPR = provider,
                    "GENCODE" = "Gencode",
                    provider
                ),
                level
            ),
            strict = FALSE
        )
        object <- new(Class = class, object)
    }
    object
}
