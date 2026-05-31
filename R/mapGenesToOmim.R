#' Map gene identifiers to OMIM disease identifiers
#'
#' Returns a named character vector (or list when multiple OMIM IDs exist per
#' gene) mapping Ensembl gene IDs or HGNC gene names to OMIM identifiers,
#' using HGNC curated data.
#'
#' @export
#' @note Updated 2026-05-31.
#'
#' @param genes `character`.
#' Human gene names (symbols) or Ensembl gene IDs.
#'
#' @param hgnc `Hgnc` or `NULL`.
#' If `NULL`, HGNC annotations will be downloaded automatically.
#'
#' @return Named `list` of `character` vectors.
#' Each element contains OMIM identifiers for the corresponding gene.
#' Genes with no OMIM mapping return `NA_character_`.
#'
#' @examples
#' x <- mapGenesToOmim(genes = c("NFE2L2", "BRCA1", "TP53"))
#' print(x)
mapGenesToOmim <-
    function(genes, hgnc = NULL) {
        if (is.null(hgnc)) {
            hgnc <- Hgnc()
        }
        assert(
            isCharacter(genes),
            is(hgnc, "Hgnc"),
            validObject(hgnc),
            isSubset(
                x = c("aliasSymbol", "ensemblGeneId", "geneName", "omimId"),
                y = colnames(hgnc)
            )
        )
        df <- as(hgnc, "DFrame")
        ## Determine whether genes look like Ensembl IDs or gene names.
        isEnsembl <- all(grepl("^ENSG[0-9]", genes))
        if (isTRUE(isEnsembl)) {
            lookupCol <- "ensemblGeneId"
            idx <- match(x = genes, table = df[[lookupCol]])
        } else {
            lookupCols <- c("geneName", "aliasSymbol", "prevSymbol")
            lookupCols <- intersect(lookupCols, colnames(df))
            idx <- matchNested(x = genes, table = df[, lookupCols])
        }
        omimIds <- df[["omimId"]]
        out <- lapply(
            X = idx,
            FUN = function(i) {
                if (is.na(i)) {
                    return(NA_character_)
                }
                ids <- omimIds[[i]]
                if (is.null(ids) || length(ids) == 0L) {
                    return(NA_character_)
                }
                as.character(ids)
            }
        )
        names(out) <- genes
        out
    }
