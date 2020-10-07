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
        cli_alert_warning(
            "{.var GRanges} does not contain biotype in {.var mcols}.\n",
            "Returning without broad class definitions."
        )
        ## Early `NA` return works successfully in `mcols()` construction.
        return(NA_character_)
    }

    ## Gene names --------------------------------------------------------------
    if ("geneName" %in% colnames(data)) {
        geneNameCol <- "geneName"
        geneNameData <- data[[geneNameCol]]
    } else {
        cli_alert_warning(
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
        cli_alert_warning("{.var GRanges} does not contain {.fun seqnames}.")
        seqnamesCol <- NULL
        seqnamesData <- NA_character_
    }

    ## Apply broad class -------------------------------------------------------
    cli_alert(sprintf(
        "Defining {.var broadClass} using: {.var %s}.",
        ## Note that `c()` call here effectively removes `NULL` definitions.
        toString(sort(c(biotypeCol, geneNameCol, seqnamesCol)))
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
