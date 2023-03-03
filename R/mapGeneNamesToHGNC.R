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
    idx <- matchNested(
        x = genes,
        table = df[, c("symbol", "prevSymbol", "aliasSymbol")]
    )
    assert(!anyNA(idx), msg = "Failed to map all genes.")
    out <- hgnc[idx, "hgncId", drop = TRUE]
    out
}
