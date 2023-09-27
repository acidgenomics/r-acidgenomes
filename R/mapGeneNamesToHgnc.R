#' Map gene names (symbols) to HGNC identifiers
#'
#' @export
#' @note Updated 2023-09-27.
#'
#' @param genes `character`.
#' Human gene names (e.g. `"TUT4"`).
#'
#' @param hgnc `Hgnc` or `NULL`.
#' If `NULL`, HGNC annotations will be downloaded automatically.
#'
#' @examples
#' ## Homo sapiens (only).
#' x <- mapGeneNamesToHgnc(genes = c("TUT4", "ZCCHC11", "TENT3A"))
#' print(x)
mapGeneNamesToHgnc <- function(genes, hgnc = NULL) {
    if (is.null(hgnc)) {
        hgnc <- Hgnc()
    }
    assert(
        isCharacter(genes),
        is(hgnc, "Hgnc")
    )
    table <- as(hgnc, "DFrame")
    table <- table[, c("symbol", "prevSymbol", "aliasSymbol")]
    idx <- matchNested(x = genes, table = table)
    assert(!anyNA(idx), msg = "Failed to map all genes.")
    out <- hgnc[idx, "hgncId", drop = TRUE]
    out
}
