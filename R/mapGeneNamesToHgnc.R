#' Map gene names (symbols) to HGNC identifiers
#'
#' @export
#' @note Updated 2025-04-14.
#'
#' @inheritParams params
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
mapGeneNamesToHgnc <-
    function(genes, ignoreCase = FALSE, hgnc = NULL) {
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
            ),
            isFlag(ignoreCase)
        )
        table <- as(hgnc, "DFrame")
        cols <- c("geneName", "prevSymbol", "aliasSymbol")
        table <- table[, cols]
        if (isTRUE(ignoreCase)) {
            genes <- toupper(genes)
            for (col in cols) {
                table[[col]] <- toupper(table[[col]])
            }
        }
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
