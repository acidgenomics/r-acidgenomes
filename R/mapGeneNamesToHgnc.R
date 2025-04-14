#' Map gene names (symbols) to HGNC identifiers
#'
#' @export
#' @note Updated 2025-04-14.
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
        is(hgnc, "Hgnc"),
        validObject(hgnc),
        isSubset(
            x = c("aliasSymbol", "geneName", "hgncId", "prevSymbol"),
            y = colnames(hgnc)
        )
    )
    table <- as(hgnc, "DFrame")
    ## Enable case-insensitive matching.
    table[["geneName2"]] <- toupper(table[["geneName"]])
    table[["prevSymbol2"]] <- toupper(table[["prevSymbol"]])
    table[["aliasSymbol2"]] <- toupper(table[["aliasSymbol"]])
    cols <- c(
        "geneName",
        "prevSymbol",
        "aliasSymbol",
        "geneName2",
        "prevSymbol2",
        "aliasSymbol2"
    )
    table <- table[, cols]
    idx <- matchNested(x = genes, table = table)
    if (anyNA(idx)) {
        fail <- genes[is.na(idx)]
        abort(sprintf(
            "%d mapping %s: %s.",
            length(fail),
            ngettext(
                n = length(fail),
                msg1 = "failure",
                msg2 = "failures"
            ),
            toInlineString(x = fail, n = 20L)
        ))
    }
    out <- hgnc[idx, "hgncId", drop = TRUE]
    out
}
