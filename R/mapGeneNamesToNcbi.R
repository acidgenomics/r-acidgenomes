#' Map gene names to NCBI
#'
#' @export
#' @note Updated 2025-04-14.
#'
#' @inheritParams AcidRoxygen::params
#' @inheritParams NcbiGeneInfo
#' @inheritParams params
#'
#' @param genes `character`.
#' Gene names (e.g. `"TUT4"`).
#'
#' @param ncbi `NcbiGeneIfo` or `NULL`.
#' If `NULL`, NCBI annotations will be downloaded automatically.
#'
#' @examples
#' ## Homo sapiens.
#' x <- mapGeneNamesToNcbi(
#'     genes = c("TUT4", "ZCCHC11", "TENT3A"),
#'     organism = "Homo sapiens"
#' )
#' print(x)
#'
#' ## Mus musculus.
#' x <- mapGeneNamesToNcbi(
#'     genes = c("Nfe2l2", "Nrf2"),
#'     organism = "Mus musculus"
#' )
#' print(x)
mapGeneNamesToNcbi <-
    function(
        genes,
        organism,
        taxonomicGroup = NULL,
        ignoreCase = FALSE,
        ncbi = NULL
    ) {
        if (is.null(ncbi)) {
            ncbi <- NcbiGeneInfo(
                organism = organism,
                taxonomicGroup = taxonomicGroup
            )
        }
        assert(
            isCharacter(genes),
            isFlag(ignoreCase),
            is(ncbi, "NcbiGeneInfo"),
            identical(organism, metadata(ncbi)[["organism"]])
        )
        table <- as(ncbi, "DFrame")
        cols <- c("geneName", "geneSynonyms")
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
        out <- ncbi[idx, "geneId", drop = TRUE]
        out <- decode(out)
        out
    }
