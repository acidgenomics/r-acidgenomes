#' Map gene names to NCBI
#'
#' @export
#' @note Updated 2025-04-07.
#'
#' @inheritParams AcidRoxygen::params
#' @inheritParams NcbiGeneInfo
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
    function(genes,
             organism,
             taxonomicGroup = NULL,
             ncbi = NULL) {
        if (is.null(ncbi)) {
            ncbi <- NcbiGeneInfo(
                organism = organism,
                taxonomicGroup = taxonomicGroup
            )
        }
        assert(
            isCharacter(genes),
            is(ncbi, "NcbiGeneInfo"),
            identical(organism, metadata(ncbi)[["organism"]])
        )
        table <- as(ncbi, "DFrame")
        table <- table[, c("geneName", "geneSynonyms")]
        idx <- matchNested(x = genes, table = table)
        if (anyNA(idx)) {
            abort(sprintf(
                "Mapping failure: %s.",
                toInlineString(genes[is.na(idx)])
            ))
        }
        out <- ncbi[idx, "geneId", drop = TRUE]
        out <- decode(out)
        out
    }
