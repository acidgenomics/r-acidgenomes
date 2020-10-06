## These code checks may be removed in a future release.
## nocov start



## Check GRanges from GFF loaded with rtracklayer against GenomicFeatures TxDb.
## Consider adding support for exons, CDS, in a future update.
.checkGRangesAgainstTxDb <- function(gr, txdb) {
    assert(is(gr, "GRanges"), is(txdb, "TxDb"))
    level <- match.arg(
        arg = metadata(gr)[["level"]],
        choices = c("genes", "transcripts")
    )
    fun <- get(
        x = level,
        envir = asNamespace("GenomicFeatures"),
        inherits = FALSE
    )
    assert(is.function(fun))
    cli_alert(sprintf("Checking %s in {.var TxDb}.", level))
    gr1 <- gr
    rm(gr)
    ## Convert the TxDb to GRanges using either `genes()` or `transcripts()`.
    ## Note that GenomicFeatures currently returns with "tx" instead of
    ## "transcript" for transcript-level annotations.
    gr2 <- fun(txdb)
    assert(is(gr2, "GRanges"))
    ## Ensure both ranges are camel case formatted.
    colnames(mcols(gr1)) <- camelCase(colnames(mcols(gr1)))
    colnames(mcols(gr2)) <- camelCase(colnames(mcols(gr2)))
    ## Always set the names on GRanges from GFF.
    idCol <- .detectGRangesIDs(gr1)
    names(gr1) <- mcols(gr1)[[idCol]]
    if (
        level == "genes" &&
        hasLength(intersect(
            x = mcols(gr1)[["geneName"]],
            y = mcols(gr2)[["geneID"]]
        ))
    ) {
        ## GenomicFeatures currently returns GFF3 input with gene symbols as the
        ## names, so ensure we're setting the GRanges from GFF to match.
        cli_alert(
            "{.var TxDb} returns gene names as identifiers for GFF3.\n",
            "Setting names on {.var GRanges} from GFF to match."
        )
        names(gr1) <- mcols(gr1)[["geneName"]]
    } else if (
        level == "transcripts" &&
        hasLength(intersect(
            x = names(gr1),
            y = mcols(gr2)[["txName"]]
        ))
    ) {
        ## `GenomicFeatures::transcripts()` returns numbers instead of correct
        ## transcript IDs, so fix that before checks. Note that this return maps
        ## the correct identifiers to `txName` column in `mcols()`.
        names(gr2) <- mcols(gr2)[["txName"]]
    } else if (
        level == "transcripts" &&
        !hasLength(intersect(
            x = names(gr1),
            y = mcols(gr2)[["txName"]]
        )) &&
        hasLength(intersect(
            x = mcols(gr1)[["transcriptName"]],
            y = mcols(gr2)[["txName"]]
        ))
    ) {
        ## GenomicFeatures currently returns GFF3 input with transcript names
        ## only, so ensure we're setting both GRanges accordingly.
        cli_alert(
            "{.var TxDb} returns gene names as identifiers for GFF3.\n",
            "Setting names on {.var GRanges} from GFF to match."
        )
        names(gr1) <- mcols(gr1)[["transcriptName"]]
        names(gr2) <- mcols(gr2)[["txName"]]
    }
    ## Ensure both GRanges are sorted by names.
    assert(hasNames(gr1), hasNames(gr2))
    gr1 <- gr1[sort(names(gr1))]
    gr2 <- gr2[sort(names(gr2))]

    ## Length and names --------------------------------------------------------
    ## Warn if GenomicFeatures has dropped ranges from GFF.
    ## This can happen for some GFF3 files and FlyBase GTF.
    if (length(gr1) > length(gr2)) {
        assert(isSubset(names(gr2), names(gr1)))
        n <- length(gr1) - length(gr2)
        cli_alert_warning(sprintf(
            fmt = paste0(
                "{.pkg GenomicFeatures} dropped %d %s from file to make ",
                "{.var TxDb}. Missing: %s."
            ),
            n,
            ngettext(
                n = n,
                msg1 = "identifier",
                msg2 = "identifiers"
            ),
            toString(setdiff(names(gr1), names(gr2)), width = 100L)
        ))
        ## Subset the GRanges from GFF to match TxDb for additional checks.
        gr1 <- gr1[names(gr2)]
    } else if (length(gr2) > length(gr1)) {
        ## This scenario happens for GENCODE GTF, where we've dropped PAR Y
        ## dupes from the GFF file but they still remain in GenomicFeatures.
        assert(isSubset(names(gr1), names(gr2)))
        n <- length(gr2) - length(gr1)
        cli_alert_warning(sprintf(
            fmt = paste0(
                "{.pkg %s} dropped %d %s from file to make {.var GRanges}. ",
                "Missing: %s."
            ),
            packageName(),
            n,
            ngettext(
                n = n,
                msg1 = "identifier",
                msg2 = "identifiers"
            ),
            toString(setdiff(names(gr1), names(gr2)), width = 100L)
        ))
        ## Subset the GRanges from TxDb to match GFF for additional checks.
        gr2 <- gr2[names(gr1)]
    }
    assert(
        identical(length(gr1), length(gr2)),
        identical(names(gr1), names(gr2))
    )

    ## Ranges and seqnames -----------------------------------------------------
    ## Compare the ranges and inform the user about mismatches.
    ## nolint start
    ## > help(topic = "IPosRanges-comparison", package = "IRanges")
    ## > vignette(topic = "IRangesOverview", package = "IRanges")
    ## nolint end
    r1 <- ranges(gr1)
    r2 <- ranges(gr2)
    diff <- r1 != r2
    if (any(diff)) {
        which <- head(which(diff), n = 10L)
        cli_alert_warning(sprintf(
            fmt = paste(
                "%d range mismatches detected in {.var TxDb}.",
                "Showing {.var GRanges} mismatch comparison (first 10).",
                "",
                "(1) {.pkg AcidGenomes}:",
                "%s",
                "",
                "(2) {.pkg GenomicFeatures}:",
                "%s",
                "",
                "If the ranges in (1) are incorrect, please file an issue:",
                "  https://github.com/acidgenomics/AcidGenomes/issues",
                "If the ranges in (2) are incorrect, please file an issue:",
                "  https://github.com/Bioconductor/GenomicFeatures/issues",
                sep = "\n"
            ),
            sum(diff, na.rm = TRUE),
            printString(r1[which]),
            printString(r2[which])
        ))
    }
    if (!identical(seqnames(gr1), seqnames(gr2))) {
        warning("'seqnames()' mismatch detected.")
    }
    invisible(TRUE)
}



#' Make TxDb from GRanges
#'
#' Note that this uses different functions depending on whether the GRanges was
#' created from a GFF3 or GTF file.
#'
#' We're now using the TxDb object for GRanges strict mode sanity checks. TxDb
#' doesn't return enough useful metadata from the original GFF/GTF file, so
#' we're parsing the return manually here using custom code instead. Note that
#' this step must be called before we attempt to sanitize any metadata columns
#' in `mcols()`.
#'
#' Note that Ensembl GFF3 currently fails this check.
#'
#' File an issue to request a fix for this in GenomicFeatures:
#'
#' `makeTxDbFromGRanges()`:
#'
#' ```
#' some exons are linked to transcripts not found in the file
#' ```
#'
#' `makeTxDbFromGFF()` works but warns:
#'
#' ```
#' The following orphan exon were dropped
#' The following orphan CDS were dropped
#' ```
#'
#' Note that `makeTxDbFromGRanges()` frequently returns warnings.
#'
#' Specifically, we're suppressing this expected warning:
#'
#' ```
#' The "phase" metadata column contains non-NA values for
#' features of type stop_codon. This information was ignored.
#' ```
#'
#' @seealso
#' - https://stackoverflow.com/questions/38603668
#' - https://stackoverflow.com/questions/16517795
#'
#' @noRd
.makeTxDbFromGFF <- function(object) {
    requireNamespaces("GenomicFeatures")
    assert(
        is(object, "GRanges"),
        isSubset("detect", names(metadata(object)))
    )
    ## Get stashed metadata values.
    file <- metadata(object)[["file"]]
    type <- metadata(object)[["detect"]][["type"]]
    txdb <- withCallingHandlers(expr = {
        ## Using `tryCatch()` here to change error message, if necessary.
        tryCatch(
            expr = if (type == "GFF3") {
                ## `makeTxDbFromGRanges()` often chokes on GRanges from GFF3,
                ## imported via `rtracklayer::import()`, so switch to using
                ## `makeTxDbFromGFF()` instead, which always works.
                cli_alert(
                    "Making {.var TxDb} using {.fun makeTxDbFromGFF}."
                )
                GenomicFeatures::makeTxDbFromGFF(file)
            } else if (type == "GTF") {
                cli_alert(
                    "Making {.var TxDb} using {.fun makeTxDbFromGRanges}."
                )
                GenomicFeatures::makeTxDbFromGRanges(object)
            },
            error = function(e) {
                stop(
                    "Failed to make TxDb using ",
                    "'GenomicFeatures::makeTxDbFromGRanges()'.\n",
                    "Set 'strict = FALSE' to disable TxDb checks.\n"
                )
            }
        )
    }, warning = function(w) {
        if (grepl(pattern = "stop_codon", x = conditionMessage(w))) {
            invokeRestart("muffleWarning")
        }
    })
    assert(is(txdb, "TxDb"))
    txdb
}



## nocov end
