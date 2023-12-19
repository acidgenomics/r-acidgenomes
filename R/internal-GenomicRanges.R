#' Add broad class annotations
#'
#' @note Updated 2023-12-19.
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
            x = names(mcols(object)),
            pattern = "^transcript"
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
    ## coercion will define `seqnames` column from the GRanges (see above).
    if ("seqnames" %in% names(df)) {
        seqnamesCol <- "seqnames"
        seqnamesData <- df[[seqnamesCol]]
    } else {
        seqnamesCol <- NULL
        seqnamesData <- NA_character_
    }
    ## Apply broad class. Note that this method doesn't seem to work right with
    ## DFrame class.
    df <- data.frame(
        "biotype" = biotypeData,
        "chromosome" = seqnamesData,
        "geneName" = geneNameData,
        row.names = names(object),
        stringsAsFactors = TRUE
    )
    x <- apply(X = df, MARGIN = 1L, FUN = .applyBroadClass)
    x <- as.factor(x)
    if (all(x == "other")) {
        alertWarning(sprintf(
            "Returning without {.var %s} in {.var %s}.",
            "broadClass", "mcols"
        ))
        return(object)
    }
    ## Disabling Rle encoding in 0.7.3. update.
    ## > x <- Rle(x)
    mcols(object)[["broadClass"]] <- x
    object
}



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
#' - *C. elegans*: Can't match by gene name. Match by `"MtDNA"` chromosome.
#' Alternatively, can match using `"MTCE"` sequence name (parent clone).
#' https://www.wormbase.org/species/c_elegans/clone/MTCE
#'
#' Note that this might not be perfect for other genomes, so consider
#' atttempting to improve support here in a future update.
#'
#' @section Ensembl biotypes:
#'
#' See [biotypes guide][] for details.
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
#' [biotypes guide]: https://m.ensembl.org/info/genome/genebuild/biotypes.html
#'
#' @section `"other"` return:
#'
#' Can return `NA_character_` here instead. Keeping this as `"other"`, to
#' main consistency with previous data sets. Also note that `NA` can behave
#' inconsistently in plotting engines.
#'
#' @seealso
#' - Can use `dplyr::case_when()` instead, which allows for a rowwise vectorized
#' if/else call stack.
#'
#' @author Rory Kirchner, Michael Steinbaugh
#' @note Updated 2023-12-19.
#' @noRd
#'
#' @param x `list`.
#' List returned via apply call using `MARGIN = 1`.
#'
#' @return `character(1)`.
.applyBroadClass <- function(x) {
    if (
        isTRUE(grepl(
            pattern = "^MT",
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
        ## Immunoglobulin.
        "ig"
    } else if (
        isTRUE(grepl(
            pattern = "^tr_",
            x = x[["biotype"]],
            ignore.case = TRUE
        ))
    ) {
        ## T cell receptor.
        "tcr"
    } else {
        "other"
    }
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



#' Include identifier version in primary identifier
#'
#' @note Updated 2021-01-30.
#' @noRd
.includeVersion <-
    function(object,
             idCol,
             idVersionCol,
             idNoVersionCol,
             quiet = TRUE) {
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



#' Match the identifier column in `GRanges` to use for names.
#'
#' @note Updated 2023-04-26.
#' @noRd
.matchGRangesNamesColumn <- function(object) {
    assert(
        is(object, "GRanges"),
        isString(metadata(object)[["level"]])
    )
    level <- match.arg(
        arg = metadata(object)[["level"]],
        choices = c("cds", "exons", "genes", "transcripts")
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



## FIXME dbxref needs to return as CompressedCharacterList.

#' Finalize `GRanges` mcols return
#'
#' @note Updated 2023-12-19.
#' @noRd
#'
#' @details
#' This step sanitizes NA values, relevels factors, and trims invalid ranges.
.returnMcols <- function(object) {
    assert(is(object, "GRanges"))
    length <- length(object)
    object <- trim(object)
    assert(hasLength(object, n = length))
    mcols <- mcols(object)
    if (isSubset("artifactualDuplication", colnames(mcols))) {
        ## e.g. GENCODE GFF.
        mcols[["artifactualDuplication"]] <- sub(
            pattern = "^real_copy_is_",
            replacement = "",
            x = mcols[["artifactualDuplication"]]
        )
    }
    if (is.character(mcols[["dbXref"]])) {
        mcols[["dbXref"]] <- CharacterList(as.list(mcols[["dbXref"]]))
    }
    if (is.list(mcols[["geneSynonyms"]])) {
        mcols[["geneSynonyms"]] <- CharacterList(mcols[["geneSynonyms"]])
    }
    if (isSubset("hgncId", colnames(mcols))) {
        mcols[["hgncId"]] <- sub(
            pattern = "^HGNC\\:",
            replacement = "",
            x = mcols[["hgncId"]]
        )
        mcols[["hgncId"]] <- as.integer(mcols[["hgncId"]])
    }
    if (is.character(mcols[["level"]])) {
        assert(allAreMatchingRegex(
            x = na.omit(mcols[["level"]]),
            pattern = "^[0-9]+$"
        ))
        mcols[["level"]] <- as.integer(mcols[["level"]])
    }
    if (is.list(mcols[["ncbiGeneId"]])) {
        mcols[["ncbiGeneId"]] <- IntegerList(mcols[["ncbiGeneId"]])
    }
    if (is.character(mcols[["partial"]])) {
        mcols[["partial"]][is.na(mcols[["partial"]])] <- "false"
        assert(isSubset(
            x = unique(mcols[["partial"]]),
            y = c("false", "true")
        ))
        mcols[["partial"]] <- as.logical(mcols[["partial"]])
    }
    if (is.character(mcols[["pseudo"]])) {
        mcols[["pseudo"]][is.na(mcols[["pseudo"]])] <- "false"
        assert(isSubset(
            x = unique(mcols[["pseudo"]]),
            y = c("false", "true")
        ))
        mcols[["pseudo"]] <- as.logical(mcols[["pseudo"]])
    }
    if (isSubset("source", colnames(mcols))) {
        ## e.g. Standardize "BestRefSeq%2CGnomon" to "BestRefSeq/Gnomon".
        mcols[["source"]] <- gsub(
            pattern = "%2C",
            replacement = "/",
            x = mcols[["source"]],
            fixed = TRUE
        )
    }
    if (
        is.character(mcols[["tag"]]) &&
        !all(is.na(mcols[["tag"]]))
    ) {
        mcols[["tag"]] <- CharacterList(as.list(mcols[["tag"]]))
    }
    if (is.integer(mcols[["txIsCanonical"]])) {
        assert(isSubset(
            x = unique(na.omit(mcols[["txIsCanonical"]])),
            y = c(0L, 1L)
        ))
        mcols[["txIsCanonical"]] <- as.logical(mcols[["txIsCanonical"]])
    }
    if (isSubset("txSupportLevel", colnames(mcols))) {
        ## Sanitize "NA (assigned to previous version 9)" to "NA".
        mcols[["txSupportLevel"]] <- sub(
            pattern = "^NA\\s.+$",
            replacement = "NA",
            x = mcols[["txSupportLevel"]]
        )
    }
    factorCols <- c(
        "exception",
        "geneBiotype",
        "geneSource",
        "level",
        "logicName",
        "seqCoordSystem",
        "source",
        "txBiotype",
        "txSource",
        "txSupportLevel",
        "type"
    )
    for (factorCol in factorCols) {
        if (
            isSubset(factorCol, colnames(mcols)) &&
            !is.factor(mcols[[factorCol]])
        ) {
            assert(is.atomic(mcols[[factorCol]]))
            mcols[[factorCol]] <- as.factor(mcols[[factorCol]])
        }
    }
    mcolsList <- lapply(
        X = mcols,
        FUN = function(x) {
            if (isS4(x) || is(x, "AsIs") || !is.atomic(x)) {
                return(x)
            }
            x <- sanitizeNa(x)
            if (all(is.na(x))) {
                return(NULL)
            }
            if (is.factor(x)) {
                x <- droplevels(x)
            }
            ## Disabling Rle encoding in 0.7.3. update.
            ## > x <- Rle(x)
            x
        }
    )
    mcolsList <- Filter(f = Negate(is.null), x = mcolsList)
    mcols <- as.DataFrame(mcolsList)
    mcols <- mcols[, sort(names(mcols)), drop = FALSE]
    mcols(object) <- mcols
    object
}



#' Standardize the `GRanges` mcols naming conventions
#'
#' @details
#' Always return using camel case, even though GFF/GTF files use snake.
#'
#' Note that this step makes `GRanges` imported via `rtracklayer::import()`
#' incompatible with `GenomicFeatures::makeTxDbFromGRanges()` parser, so be
#' sure to call that function prior to attempting to run this step.
#'
#' @note Updated 2023-12-19.
#' @noRd
.standardizeMcolsNames <- function(object) {
    assert(is(object, "GRanges"))
    mcols <- mcols(object)
    names(mcols) <- camelCase(names(mcols), strict = TRUE)
    ## Ensure "tx" prefix is used consistently instead of "transcript".
    ## This convention was changed in v0.2.0 release.
    names(mcols) <- gsub(
        pattern = "^transcript",
        replacement = "tx",
        x = names(mcols),
        ignore.case = FALSE
    )
    ## Ensure "Id" is always capitalized at the end.
    names(mcols) <-
        gsub(
            pattern = "(.+)id$",
            replacement = "\\1Id",
            x = names(mcols),
            ignore.case = FALSE
        )
    ## Use NCBI instead of Entrez for database name.
    if (isSubset("entrezId", names(mcols))) {
        names(mcols)[names(mcols) == "entrezId"] <- "ncbiGeneId"
    }
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
    ## e.g. GENCODE GFF.
    names(mcols)[names(mcols) == "artifDupl"] <- "artifactualDuplication"
    names(mcols)[names(mcols) == "dbxref"] <- "dbXref"
    mcols(object) <- mcols
    object
}



#' Make genomic ranges (`GRanges`)
#'
#' This is the main `GRanges` final return generator, used by
#' `makeGRangesFromEnsembl()` and `makeGRangesFromGff()`.
#'
#' @note Updated 2023-12-19.
#' @noRd
.makeGRanges <-
    function(object,
             ignoreVersion,
             extraMcols) {
        assert(
            is(object, "GRanges"),
            hasLength(object),
            isFlag(ignoreVersion),
            isFlag(extraMcols),
            isString(metadata(object)[["level"]]),
            isString(metadata(object)[["provider"]])
        )
        level <- match.arg(
            arg = metadata(object)[["level"]],
            choices = .grangesLevels
        )
        provider <- metadata(object)[["provider"]]
        object <- .standardizeMcolsNames(object)
        if (isFALSE(ignoreVersion)) {
            object <- .includeGeneVersion(object)
            object <- .includeTxVersion(object)
        }
        if (isTRUE(extraMcols)) {
            object <- .addBroadClass(object)
            if (isSubset(
                x = metadata(object)[["provider"]],
                y = c("Ensembl", "GENCODE")
            )) {
                object <- .addEnsemblFtpMcols(
                    object = object,
                    ignoreVersion = ignoreVersion
                )
            }
        }
        object <- .returnMcols(object)
        idCol <- .matchGRangesNamesColumn(object)
        assert(isSubset(idCol, names(mcols(object))))
        alert(sprintf(
            "Defining names by {.var %s} column in {.fun %s}.",
            idCol, "mcols"
        ))
        ## Sort the ranges by genomic location.
        object <- sort(object)
        ## Sort the metadata elements alphabetically.
        metadata(object) <- append(
            x = metadata(object),
            values = list(
                "date" = Sys.Date(),
                "ignoreVersion" = ignoreVersion,
                "packageVersion" = .pkgVersion
            )
        )
        metadata(object) <- metadata(object)[sort(names(metadata(object)))]
        if (identical(provider, "RefSeq")) {
            alertInfo(sprintf(
                "Splitting {.cls %s} by {.var %s} into {.cls %s}.",
                "GRanges", idCol, "GRangesList"
            ))
            ## Metadata gets dropped during `split()` call; stash and reassign.
            meta <- metadata(object)
            object <- split(x = object, f = as.factor(mcols(object)[[idCol]]))
            metadata(object) <- meta
        } else {
            names <- as.character(mcols(object)[[idCol]])
            assert(
                hasNoDuplicates(names),
                !anyNA(names),
                msg = paste(
                    "Invalid and/or duplicated identifiers detected.",
                    "Setting 'ignoreVersion = FALSE' may resolve this."
                )
            )
            names(object) <- names
            ## This check fails for split GRangesList.
            assert(
                isFALSE(is.unsorted(object)),
                msg = "Ranges are not sorted."
            )
        }
        ## Run final assert checks before returning.
        assert(validObject(object))
        if (isSubset(level, c("genes", "transcripts"))) {
            class <- upperCamelCase(
                object = tolower(paste(provider, level)),
                strict = FALSE
            )
            object <- new(Class = class, object)
        }
        object
    }
