#' Map gene names to NCBI
#'
#' @export
#' @note Updated 2023-03-03.
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
#' x <- mapGeneNamesToNCBI(
#'     genes = c("TUT4", "ZCCHC11", "TENT3A"),
#'     organism = "Homo sapiens"
#' )
#' print(x)
#'
#' ## Mus musculus.
#' x <- mapGeneNamesToNCBI(
#'     genes = c("Nfe2l2", "Nrf2"),
#'     organism = "Mus musculus"
#' )
#' print(x)
mapGeneNamesToNCBI <-
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
        table <- as(ncbi, "DataFrame")
        table <- table[, c("geneName", "geneSynonyms")]
        idx <- matchNested(x = genes, table = table)
        assert(!anyNA(idx), msg = "Failed to map all genes.")
        out <- ncbi[idx, "geneId", drop = TRUE]
        out <- decode(out)
        out
    }
