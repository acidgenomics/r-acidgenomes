#' Map human gene names to mouse orthologs
#'
#' Useful for PDX (patient-derived xenograft) databases where human gene
#' expression needs to be mapped to mouse reference annotations.
#'
#' @export
#' @note Updated 2026-05-31.
#'
#' @param genes `character`.
#' Human gene names (symbols) to map to mouse orthologs.
#'
#' @param jax `JaxHumanToMouse` or `NULL`.
#' If `NULL`, JAX human-to-mouse ortholog data will be downloaded
#' automatically.
#'
#' @param column `character(1)`.
#' Which mouse identifier column to return.
#' One of `"mouseGeneName"`, `"mouseNcbiGeneId"`, or `"mouseMgiId"`.
#'
#' @return Named `character` (or `integer` for NCBI/MGI columns).
#' Mouse ortholog identifiers, named by the input human gene names.
#' Genes with no ortholog return `NA`.
#'
#' @examples
#' x <- mapHumanToMouse(genes = c("TP53", "BRCA1", "NFE2L2"))
#' print(x)
mapHumanToMouse <-
    function(
        genes,
        jax = NULL,
        column = c("mouseGeneName", "mouseNcbiGeneId", "mouseMgiId")
    ) {
        if (is.null(jax)) {
            jax <- JaxHumanToMouse()
        }
        assert(
            isCharacter(genes),
            is(jax, "JaxHumanToMouse"),
            validObject(jax)
        )
        column <- match.arg(column)
        df <- as(jax, "DFrame")
        assert(isSubset(
            x = c("humanGeneName", column),
            y = colnames(df)
        ))
        idx <- match(x = genes, table = df[["humanGeneName"]])
        out <- df[[column]][idx]
        names(out) <- genes
        out
    }
