#' Update gene symbols
#'
#' @export
#' @note Updated 2023-12-21.
#'
#' @param geneNames `character`.
#' Gene names (symbols).
#'
#' @param organism `character(1)`.
#' Full Latin organism name.
#' Human genes are checked against HGNC.
#' Other organisms are checked against NCBI.
#'
#' @return `character`.
#' Named character, with original input defined in the names.
#'
#' @seealso
#' - https://twitter.com/samuel_marsh/status/1737855289876201922
#' - https://waldronlab.io/HGNChelper/
#'
#' @examples
#' ## Homo sapiens ====
#' object <- updateGeneSymbols(
#'     geneNames = c("ZCCHC6", "ZCCHC11"),
#'     organism = "Homo sapiens"
#' )
#' print(object)
#'
#' ## Mus musculus ====
#' object <- updateGeneSymbols(
#'     geneNames = c("Zcchc6", "Zcchc11"),
#'     organism = "Mus musculus"
#' )
#' print(object)
updateGeneSymbols <- function(geneNames, organism) {
    assert(
        isCharacter(geneNames),
        hasNoDuplicates(geneNames),
        !anyNA(geneNames),
        isOrganism(organism)
    )
    switch(
        EXPR = organism,
        "Homo sapiens" = {
            df <- Hgnc()
            cols <- c("geneName", "aliasSymbol", "prevSymbol")
        },
        {
            df <- NcbiGeneInfo(organism = organism)
            cols <- c("geneName", "geneSynonyms")
        }
    )
    df <- as(df, "DFrame")
    assert(isSubset(cols, colnames(df)))
    df <- df[, cols]
    idx <- matchNested(geneNames, df)
    if (anyNA(idx)) {
        fail <- geneNames[which(is.na(idx))]
        abort(sprintf(
            "Failed to match %s gene %s: %s.",
            length(fail),
            ngettext(
                n = length(fail),
                msg1 = "symbol",
                msg2 = "symbols"
            ),
            toInlineString(fail)
        ))
    }
    out <- df[["geneName"]][idx]
    names(out) <- geneNames
    out
}
