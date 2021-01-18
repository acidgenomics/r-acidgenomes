## GFF metadata detection ======================================================
#' Detect the GFF source information
#'
#' @details
#' Assuming we've already cached the URL using BiocFileCache here.
#' This step will load into GRanges via rtracklayer.
#'
#' @note Updated 2021-01-18.
#' @noRd
.detectGFF <- function(object) {
    assert(is(object, "GRanges"))
    c(
        "source" = .detectGFFSource(object),
        "type" = .detectGFFType(object)
    )
}



#' Detect the database source of the genome annotations
#'
#' @note Updated 2021-01-18.
#' @noRd
.detectGFFSource <- function(object) {
    assert(is(object, "GRanges"))
    mcols <- mcols(object)
    source <- mcols[["source"]]
    if (
        ## UCSC (e.g. hg38_knownGene)
        any(grepl(pattern = "_knownGene$", x = source, ignore.case = FALSE))
    ) {
        ## nocov start
        stop(paste0(
            "UCSC genomes are intentionally not supported.\n",
            "Use a pre-built TxDb package instead ",
            "(e.g. 'TxDb.Hsapiens.UCSC.hg38.knownGene')."
        ))
        ## nocov end
    } else if (
        ## Check for GENCODE prior to Ensembl.
        any(source == "ENSEMBL") &&
        any(source == "HAVANA") &&
        "gene_type" %in% colnames(mcols)
    ) {
        out <- "GENCODE"
    } else if (
        any(grepl(pattern = "FlyBase", x = source, ignore.case = FALSE))
    ) {
        out <- "FlyBase"
    } else if (
        any(grepl(pattern = "WormBase", x = source, ignore.case = FALSE))
    ) {
        out <- "WormBase"
    } else if (
        any(grepl(pattern = "RefSeq", x = source, ignore.case = FALSE))
    ) {
        out <- "RefSeq"
    } else if (
        any(grepl(
            pattern = "ensembl|havana",
            x = source,
            ignore.case = FALSE
        ))
    ) {
        out <- "Ensembl"
    } else {
        ## nocov start
        stop(sprintf(
            fmt = paste(
                "Failed to detect valid GFF/GTF source.",
                "Supported: %s",
                sep = "\n"
            ),
            toString(c("Ensembl", "FlyBase", "GENCODE", "RefSeq", "WormBase"))
        ))
        ## nocov end
    }
    out
}



#' Determine if input is GFF3 or GTF (GFF2)
#'
#' @note Updated 2021-01-18.
#' @noRd
.detectGFFType <- function(object) {
    assert(is(object, "GRanges"))
    if (any(c("ID", "Name", "Parent") %in% colnames(mcols(object)))) {
        out <- "GFF3"
    } else {
        out <- "GTF"
    }
    out
}



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
#' @note Updated 2021-01-18.
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
        ## FlyBase GTF will hit this step.
        ## Note that we're early returning without calculations in this case.
        alertWarning(paste(
            "{.var GRanges} does not contain biotype in {.var mcols}.",
            "Returning without broad class definitions."
        ))
        ## Early `NA` return works successfully in `mcols()` construction.
        return(NA_character_)
    }
    ## Gene names --------------------------------------------------------------
    if ("geneName" %in% colnames(df)) {
        geneNameCol <- "geneName"
        geneNameData <- df[[geneNameCol]]
    } else {
        alertWarning(
            "{.var GRanges} does not contain gene names in {.fun mcols}."
        )
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
        ## Don't think this is currently possible to hit, but keep just in case.
        alertWarning("{.var GRanges} does not contain {.fun seqnames}.")
        seqnamesCol <- NULL
        seqnamesData <- NA_character_
    }
    ## Apply broad class -------------------------------------------------------
    alert(sprintf(
        "Defining {.var broadClass} using: %s.",
        ## Note that `c()` call here effectively removes `NULL` definitions.
        toInlineString(sort(c(biotypeCol, geneNameCol, seqnamesCol)))
    ))
    ## Note that this method doesn't seem to work right with DataFrame class.
    df <- data.frame(
        "biotype" = biotypeData,
        "chromosome" = seqnamesData,
        "geneName" = geneNameData,
        stringsAsFactors = TRUE
    )
    ## Consider adding BiocParallel support here for improved speed.
    mcols(object)[["broadClass"]] <- Rle(as.factor(apply(
        X = df,
        MARGIN = 1L,
        FUN = .applyBroadClass
        )))
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
#' Append the transcript version to the identifier (e.g. ENST00000000233.10).
#'
#' @note Updated 2021-01-18.
#' @noRd
.addGeneVersion <- function(object) {
    assert(is(object, "GRanges"))
    if (!isSubset("geneId", colnames(mcols(object)))) {
        return(object)
    }
    alert("Adding version to gene identifiers.")
    mcolnames <- colnames(mcols(object))
    if (isSubset("geneIdVersion", mcolnames)) {
        ## `makeGRangesFromEnsembl()` output via ensembldb.
        id <- mcols(object)[["geneIdVersion"]]
    } else if (isSubset("geneVersion", mcolnames)) {
        ## FIXME NEED TO RECHECK THIS FOLLOWING TxDb approach update.
        ## `makeGRangesFromGFF()` output.
        id <- mcols(object)[["geneId"]]
        version <- mcols(object)[["geneVersion"]]
        id <- Rle(paste(id, version, sep = "."))
    } else {
        stop("Failed to locate gene identifier version.")  # nolint
    }
    mcols(object)[["geneId"]] <- id
    object
}



## FIXME This needs to be called when `ignoreVersion` is changed.
## FIXME WE SHOULD ALSO MODIFY THE GENE VERSION HERE, NO?
## FIXME THIS NEEDS TO MATCH UP WITH TX2GENE OUTPUT...

#' Add the transcript identifier version
#'
#' Append the transcript version to the identifier (e.g. ENST00000000233.10).
#'
#' @note Updated 2021-01-18.
#' @noRd
.addTxVersion <- function(object) {
    assert(
        is(object, "GRanges"),
        isSubset("txId", colnames(mcols(object)))
    )
    alert("Adding version to transcript identifiers.")
    mcolnames <- colnames(mcols(object))
    if (isSubset("txIdVersion", mcolnames)) {
        ## `makeGRangesFromEnsembl()` output via ensembldb.
        id <- mcols(object)[["txIdVersion"]]
    } else if (isSubset("txVersion", mcolnames)) {
        ## `makeGRangesFromGFF()` output.
        id <- mcols(object)[["txId"]]
        version <- mcols(object)[["txVersion"]]
        id <- Rle(paste(id, version, sep = "."))
    } else {
        stop("Failed to locate transcript identifier version.")  # nolint
    }
    mcols(object)[["txId"]] <- id
    object
}



## Source-specific metadata parsers ============================================
## FIXME ADD makeGRangesFromWormBase
## FIXME ADD makeGRangesFromFlyBase
## FIXME ADD makeGRangesFromRefSeq
##
## FIXME NEEDS TO REIMPORT GTF/GFF INTO GRANGES AND WORK FROM THERE.
## FIXME (each of these need to support level argument)

## Legacy code:
## https://github.com/acidgenomics/r-acidgenomes/blob/main/R/makeGRangesFromGFF-FlyBase.R
## https://github.com/acidgenomics/r-acidgenomes/blob/main/R/makeGRangesFromGFF-RefSeq.R
## https://github.com/acidgenomics/r-acidgenomes/blob/main/R/makeGRangesFromGFF-WormBase.R



## Standardization =============================================================
#' Match the identifier column in GRanges to use for names.
#'
#' @note Updated 2021-01-18.
#' @noRd
.matchGRangesNamesColumn <- function(object) {
    assert(is(object, "GRanges"))
    level <- match.arg(
        arg = metadata(object)[["level"]],
        choices = c("genes", "transcripts")
    )
    x <- switch(
        EXPR = level,
        "genes" = "geneId",
        "transcripts" = "txId"
    )
    assert(isSubset(x, colnames(mcols(object))))
    x
}



#' Minimize GRanges mcols
#'
#' This step sanitizes NA values, applies run-length encoding (to reduce memory
#' overhead), and trims any invalid ranges.
#'
#' @note Updated 2021-01-18.
#' @noRd
.minimizeGRanges <- function(object) {
    assert(is(object, "GRanges"))
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
                if (is.factor(x)) {
                    x <- droplevels(x)
                }
                x <- Rle(x)
                x
            }
        }
    )
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
#' @note Updated 2021-01-18.
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
        alert(
            "Renaming {.var symbol} to {.var geneName} in {.fun mcols}."
        )
        mcols[["geneName"]] <- mcols[["symbol"]]
        mcols[["symbol"]] <- NULL
    }
    mcols(object) <- mcols
    object
}



## Main generator ==============================================================
## FIXME ignoreVersion = FALSE needs to return gene ID versions.
## FIXME NEED TO IMPROVE CONSISTENCY Of METADATA RETURN.
## FIXME NEED TO DETECT ORGANISM AUTOMATICALLY.
## FIXME NEED TO SLOT ARGUMENTS INTO OBJECT.
## FIXME ORGANISM ALWAYS NEEDS TO BE DEFINED.
## FIXME DONT ALLOW GRANGESLIST RETURN HERE.
## FIXME RETURN THE OBJECTS CLASSED BY GENOME.
## FIXME HANDLE THE TXVERSION DIFFERENTLY HERE?
## FIXME SWITCH FROM "TRANSCRIPT" PREFIX TO "TX", FOR BIOC CONSISTENCY.
## FIXME Rework how we return these objects classed by source and level.

#' Make GRanges
#'
#' This is the main GRanges final return generator, used by
#' `makeGRangesFromEnsembl()` and `makeGRangesFromGFF()`.
#'
#' @note Updated 2021-01-14.
#' @noRd
.makeGRanges <- function(
    object,
    ignoreVersion = TRUE,
    broadClass = TRUE,
    synonyms = TRUE
) {
    assert(
        is(object, "GRanges"),
        hasNames(object),
        hasLength(object),
        isFlag(ignoreVersion),
        isFlag(broadClass),
        isFlag(synonyms)
    )
    object <- .minimizeGRanges(object)
    object <- .standardizeGRanges(object)
    ## FIXME RETHINK THESE STEPS <<<
    if (isFALSE(ignoreVersion)) {
        ## FIXME REWORK THIS STEP.
        object <- .addGeneVersion(object)
        object <- .addTxVersion(object)
    }
    if (isTRUE(broadClass)) {
        object <- .addBroadClass(object)
    }
    if (isTRUE(synonyms)) {
        object <- .addGeneSynonyms(object)
    }


    ## FIXME RETHINK THIS APPROACH.
    ## Ensure the ranges are sorted by identifier.
    ## FIXME NEED TO APPLY IDENTIFIER VERSION HANDLING PRIOR TO THIS STEP.
    idCol <- .matchGRangesNamesColumn(object)
    alert(sprintf("Arranging by {.var %s}.", idCol))
    names(object) <- mcols(object)[[idCol]]
    object <- object[sort(names(object))]





    ## Sort the metadata columns alphabetically.
    mcols(object) <-
        mcols(object)[, sort(colnames(mcols(object))), drop = FALSE]
    ## Prepare the metadata.
    ## Slot organism into metadata.
    ## FIXME SHOULD WE RETHINK THIS STEP?
    object <- .slotOrganism(object)
    ## Ensure object contains prototype metadata.
    metadata(object) <- c(.prototypeMetadata, metadata(object))
    idCol <- .matchGRangesNamesColumn(object)
    assert(isSubset(idCol, colnames(mcols(object))))
    names <- as.character(mcols(object)[[idCol]])
    assert(!any(is.na(names)))
    ## Inform the user if the object contains invalid names, showing offenders.
    invalid <- setdiff(names, make.names(names, unique = TRUE))
    if (hasLength(invalid)) {
        invalid <- sort(unique(invalid))
        alertWarning(sprintf(
            fmt = "%d invalid %s: %s.",
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
    ## FIXME RETHINK THIS, UNNECESSARY ONCE WE FIX REFSEQ...
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
    ## FIXME ADD ASSERT FOR ORGANISM IN METADATA HERE.
    assert(!any(is.na(seqlengths(object))))
    validObject(object)
    object
}
