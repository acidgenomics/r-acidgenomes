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
#' @note Updated 2020-01-20.
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
        ## immunoglobulin
        "ig"
    } else if (
        isTRUE(grepl(
            pattern = "^tr_",
            x = x[["biotype"]],
            ignore.case = TRUE
        ))
    ) {
        ## T cell receptor
        "tcr"
    } else {
        "other"
    }
}



#' Broad class definitions factor return
#'
#' @note Updated 2020-10-07.
#' @noRd
#'
#' @inheritParams AcidRoxygen::params
#' @param object `GRanges`.
#'
#' @return `factor`.
.broadClass <- function(object) {
    assert(is(object, "GRanges"))
    colnames(mcols(object)) <- camelCase(colnames(mcols(object)))
    ## Early return if already defined in `mcols()`.
    if ("broadClass" %in% colnames(mcols(object))) {
        out <- mcols(object)[["broadClass"]]
        out <- as.factor(out)
        return(out)
    }
    ## Need to strip the names on the object here, otherwise data.frame coercion
    ## will error if the object contains duplicate names, which can happen with
    ## GRanges that need to be split to GRangesList.
    names(object) <- NULL
    ## This step coerces the genomic ranges, including associated metadata
    ## stored in `mcols()` to a flat data frame, which will now contain
    ## seqnames, start, end, width, and strand as columns.
    data <- as.data.frame(object)
    ## Biotypes ----------------------------------------------------------------
    ## Prioritizing transcript biotype over gene, if defined. This only applies
    ## for transcript-level GRanges. For gene-level GRanges, the gene biotypes
    ## will be used, as expected.
    if ("transcriptBiotype" %in% colnames(data)) {
        biotypeCol <- "transcriptBiotype"
        biotypeData <- data[[biotypeCol]]
    } else if ("geneBiotype" %in% colnames(data)) {
        biotypeCol <- "geneBiotype"
        biotypeData <- data[[biotypeCol]]
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
    if ("geneName" %in% colnames(data)) {
        geneNameCol <- "geneName"
        geneNameData <- data[[geneNameCol]]
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
    if ("seqnames" %in% colnames(data)) {
        seqnamesCol <- "seqnames"
        seqnamesData <- data[[seqnamesCol]]
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
    data <- data.frame(
        biotype = biotypeData,
        chromosome = seqnamesData,
        geneName = geneNameData,
        stringsAsFactors = TRUE
    )
    ## Consider adding BiocParallel support here for improved speed.
    out <- apply(X = data, MARGIN = 1L, FUN = .applyBroadClass)
    out <- as.factor(out)
    out
}



#' Detect the GFF source information
#'
#' @details
#' Assuming we've already cached the URL using BiocFileCache here.
#' This step will load into GRanges via rtracklayer.
#'
#' @note Updated 2021-01-12.
#' @noRd
.detectGFF <- function(object) {
    assert(is(object, "GRanges"))
    alert("Detecting annotation source.")
    source <- .detectGFFSource(object)
    type <- .detectGFFType(object)
    alertInfo(sprintf("%s %s detected.", source, type))
    out <- c(
        "source" = source,
        "type" = type
    )
    out
}



#' Report the source of the gene annotations
#'
#' @note Updated 2021-01-10.
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



#' Determine if GFF or GTF
#'
#' @note Updated 2021-01-10.
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



## FIXME ignoreVersion = FALSE needs to return gene ID versions.
## FIXME NEED TO IMPROVE CONSISTENCY Of METADATA RETURN.
## FIXME NEED TO DETECT ORGANISM AUTOMATICALLY.
## FIXME NEED TO SLOT ARGUMENTS INTO OBJECT.
## FIXME ORGANISM ALWAYS NEEDS TO BE DEFINED.
## FIXME DONT ALLOW GRANGESLIST RETURN HERE.
## FIXME RETURN THE OBJECTS CLASSED BY GENOME.
## FIXME HANDLE THE TXVERSION DIFFERENTLY HERE?

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
    object <- trim(object)
    length <- length(object)
    object <- .minimizeGRanges(object)
    object <- .standardizeGRanges(object)
    assert(hasLength(object, n = length))
    if (isFALSE(ignoreVersion)) {
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
    ## FIXME SHOULD WE RETHINK THIS STEP?
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
