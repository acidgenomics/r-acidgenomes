## FIXME RETHINK METADATA RETURN STRUCTURE



#' Make GRanges from EnsDb object
#'
#' @details
#' Use specific `EnsDb` object as annotation source.
#' Alternatively, can pass in an EnsDb package name as a `character(1)`.
#'
#' @export
#' @note Updated 2020-10-06.
#'
#' @inheritParams AcidRoxygen::params
#' @inheritParams params
#'
#' @return `GRanges`.
#'
#' @examples
#' if ("EnsDb.Hsapiens.v75" %in% rownames(installed.packages())) {
#'     x <- makeGRangesFromEnsDb("EnsDb.Hsapiens.v75")
#' }
makeGRangesFromEnsDb <- function(
    object,
    level = c("genes", "transcripts"),
    ignoreTxVersion = TRUE,
    broadClass = TRUE,
    synonyms = FALSE
) {
    assert(
        isFlag(ignoreTxVersion),
        isFlag(broadClass),
        isFlag(synonyms)
    )
    level <- match.arg(level)
    alert("Making {.var GRanges} from {.var EnsDb}.")
    userAttached <- .packages()
    if (isString(object)) {
        package <- object
        requireNamespaces(package)
        object <- get(
            x = package,
            envir = asNamespace(package),
            inherits = FALSE
        )
    }
    assert(is(object, "EnsDb"))
    metadata <- .getEnsDbMetadata(object, level = level)
    genes <- genes(
        x = object,
        order.by = "gene_id",
        return.type = "GRanges"
    )
    if (level == "genes") {
        out <- genes
    }
    if (level == "transcripts") {
        ## FIXME RETHINK THIS...
        transcripts <- transcripts(
            x = object,
            order.by = "tx_id",
            return.type = "GRanges"
        )
        ## FIXME RETHINK THIS APPROACH.
        out <- .mergeGenesIntoTranscripts(
            transcripts = transcripts,
            genes = genes
        )
    }
    forceDetach(keep = userAttached)
    metadata(out) <- metadata
    out <- .makeGRanges(
        object = out,
        ignoreTxVersion = ignoreTxVersion,
        broadClass = broadClass,
        synonyms = synonyms
    )
    out
}
