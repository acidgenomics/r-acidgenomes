#' Map gene names (symbols) to HGNC identifiers
#'
#' @export
#' @note Updated 2023-03-03.
#'
#' @param genes `character`.
#' Gene names (e.g. `"TUT4"`).
#'
#' @param hgnc `HGNC` or `NULL`.
#' If `NULL`, HGNC annotations will be downloaded automatically.
#'
#' @examples
#' x <- mapGeneNamesToHGNC(genes = c("TUT4", "ZCCHC11", "TENT3A"))
mapGeneNamesToHGNC <- function(genes, hgnc = NULL) {
    if (is.null(hgnc)) {
        hgnc <- HGNC()
    }
    assert(
        isCharacter(genes),
        is(hgnc, "HGNC")
    )
    df <- as(hgnc, "DataFrame")
    ## FIXME Move this to AcidBase.
    matchNested <- function(x, table) {
        lst <- apply(
            X = table,
            MARGIN = 1L,
            FUN = function(x) {
                x <- unlist(x, recursive = TRUE, use.names = FALSE)
                x <- na.omit(x)
                x <- unique(x)
                x
            },
            simplify = FALSE
        )
        idx <- rep(
            x = seq_along(lst),
            times = vapply(
                X = lst,
                FUN = length,
                FUN.VALUE = integer(1L)
            )
        )
        value <- unlist(x = lst, recursive = FALSE, use.names = FALSE)
        df <- data.frame("idx" = idx, "value" = value)
        df <- df[!duplicated(df[["value"]]), , drop = FALSE]
        idx <- match(x = x, table = df[["value"]])
        out <- df[["idx"]][idx]
        out
    }
    idx <- matchNested(
        x = genes,
        table = df[, c("symbol", "prevSymbol", "aliasSymbol")]
    )
    assert(!anyNA(idx), msg = "Failed to map all genes.")
    out <- hgnc[idx, "hgncId", drop = TRUE]
    out
}
