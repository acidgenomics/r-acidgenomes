#' Map gene names (symbols) to HGNC identifiers
#'
#' @export
#' @note Updated 2023-03-02.
#'
#' @param genes `character`.
#' Gene names (e.g. `"TP53"`).
#'
#' @param hgnc `HGNC` or `NULL`.
#' If `NULL`, `HGNC` annotations will be downloaded automatically.
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
    hgncDf <- as(hgnc, "DataFrame")
    pool <- function(df) {
        lst <- apply(
            X = df,
            MARGIN = 1L,
            FUN = function(x) {
                x <- unlist(x, recursive = TRUE, use.names = FALSE)
                x <- na.omit(x)
                x <- unique(x)
                x
            },
            simplify = FALSE
        )
        rep <- rep(
            x = seq_along(lst),
            times = vapply(
                X = lst,
                FUN = length,
                FUN.VALUE = integer(1L)
            )
        )
        unlist <- unlist(x = lst, recursive = FALSE, use.names = FALSE)
        out <- list("rep" = rep, "unlist" = unlist)
        out
    }
    pool <- pool(df = hgncDf[, c("symbol", "prevSymbol", "aliasSymbol")])
    idx <- match(x = genes, table = pool[["unlist"]])
    assert(
        !anyNA(idx),
        msg = "Failed to map all genes."
    )
    idx2 <- pool[["rep"]][idx]
    out <- hgnc[idx2, "hgncId", drop = TRUE]
    out
}
