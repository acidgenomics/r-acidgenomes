## FIXME RETHINK METADATA RETURN STRUCTURE.
## FIXME RETHINK ALLOWING BROADCLASS AND SYNONYMS OPTIONS HERE.



#' Make GRanges from EnsDb object
#'
#' @details
#' Use specific `EnsDb` object as annotation source.
#' Alternatively, can pass in an EnsDb package name as a `character(1)`.
#'
#' @export
#' @note Updated 2021-01-10.
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
    ignoreVersion = TRUE,  # FIXME RENAME THIS.
    ## FIXME CONSIDER RENAMING THIS TO SIMPLY FAST.
    broadClass = TRUE,
    synonyms = FALSE
) {
    assert(
        isFlag(ignoreVersion),
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
    args <- list(
        "x" = object,
        "return.type" = "GRanges"
    )
    geneCols <- c(
        "gene_id",
        "gene_name",
        "gene_biotype",
        "seq_coord_system",
        "entrezid"
    )
    switch(
        EXPR = level,
        "genes" = {
            fun <- ensembldb::genes
            args <- append(
                x = args,
                values = list(
                    "columns" = geneCols,
                    "order.by" = "gene_id",
                    "single.strand.genes.only" = TRUE
                )
            )
        },
        "transcripts" = {
            fun <- ensembldb::transcripts
            args <- append(
                x = args,
                values = list(
                    "columns" = c(
                        "tx_id",
                        "tx_name",
                        "tx_biotype",
                        "tx_cds_seq_start",
                        "tx_cds_seq_end",
                        geneCols
                    ),
                    "order.by" = "tx_id"
                )
            )
        }
    )
    ## This step can warn about out-of-bound ranges that need to be trimmed.
    suppressWarnings({
        gr <- do.call(what = fun, args = args)
    })
    assert(is(gr, "GRanges"))
    forceDetach(keep = userAttached)
    metadata(gr) <- metadata
    .makeGRanges(
        object = gr,
        ## FIXME RENAME TO IGNOREVERSION?
        ignoreVersion = ignoreVersion,
        ## FIXME SIMPLIFY THE PASSTHROUGH OF THIS...
        broadClass = broadClass,
        synonyms = synonyms
    )
}
