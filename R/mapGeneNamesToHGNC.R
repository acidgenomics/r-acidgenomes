#' Map gene names (symbols) to HGNC identifiers
#'
#' @export
#' @note Updated 2023-03-03.
#'
#' @param genes `character`.
#' Human gene names (e.g. `"TUT4"`).
#'
#' @param hgnc `HGNC` or `NULL`.
#' If `NULL`, HGNC annotations will be downloaded automatically.
#'
#' @examples
#' ## Homo sapiens (only).
#' x <- mapGeneNamesToHGNC(genes = c("TUT4", "ZCCHC11", "TENT3A"))
#' print(x)
mapGeneNamesToHGNC <- function(genes, hgnc = NULL) {
    if (is.null(hgnc)) {
        hgnc <- HGNC()
    }
    assert(
        isCharacter(genes),
        is(hgnc, "HGNC")
    )
    table <- as(hgnc, "DFrame")
    table <- table[, c("symbol", "prevSymbol", "aliasSymbol")]
    idx <- matchNested(x = genes, table = table)
    assert(!anyNA(idx), msg = "Failed to map all genes.")
    out <- hgnc[idx, "hgncId", drop = TRUE]
    out
}
