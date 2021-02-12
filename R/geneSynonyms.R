#' Gene synonyms
#'
#' Look up gene synonyms from NCBI.
#'
#' @note Updated 2021-02-12.
#' @export
#'
#' @section *Caenorhabditis elegans*:
#'
#' Synonym support for *Caenorhabditis elegans* is poor on NCBI. Use the
#' [WormBase](https://r.acidgenomics.com/packages/wormbase/) package instead.
#'
#'
#' @inheritParams AcidRoxygen::params
#'
#' @return `DataFrame` containing `geneId` and `geneSynonyms` columns.
#'
#' @examples
#' object <- geneSynonyms(
#'     organism = "Homo sapiens",
#'     taxonomicGroup = "Mammalia",
#'     geneIDType = "Ensembl"
#' )
#' print(object)
geneSynonyms <- function(
    organism,
    taxonomicGroup = NULL,
    geneIDType = c("Entrez", "Ensembl", "HGNC", "OMIM")
) {
    geneIDType <- match.arg(geneIDType)
    df <- EntrezGeneInfo(
        organism = organism,
        taxonomicGroup = taxonomicGroup
    )
    cols <- c("geneId", "geneSynonyms", "dbXrefs")
    assert(
        isSubset(cols, colnames(df)),
        is(df[["geneSynonyms"]], "CharacterList")
    )
    df <- as(df, "DataFrame")
    df <- df[, cols, drop = FALSE]
    keep <- !all(is.na(df[["geneSynonyms"]]))
    df <- df[keep, , drop = FALSE]
    if (identical(geneIDType, "Entrez")) {
        df[["geneId"]] <- decode(df[["geneId"]])
        assert(is.integer(df[["geneId"]]))
        df[["dbXrefs"]] <- NULL
        return(df)
    }
    pattern <- paste0(
        "^",
        switch(
            EXPR = geneIDType,
            "HGNC" = "HGNC:HGNC",
            "OMIM" = "MIM",
            geneIDType
        ),
        ":(.+)$"
    )
    keep <- any(grepl(pattern = pattern, x = df[["dbXrefs"]]))
    df <- df[keep, , drop = FALSE]
    ## Extract the gene identifiers.
    x <- df[["dbXrefs"]]
    x <- x[grepl(pattern = pattern, x = x)]
    x <- gsub(pattern = pattern, replacement = "\\1", x = x)
    assert(is(x, "CharacterList"))
    ## Handle cases where identifiers don't map 1:1 to Entrez.
    ## In this case, the first (oldest) identifier will be used.
    x <- unlist(
        x = lapply(X = x, FUN = `[[`, 1L),
        recursive = FALSE,
        use.names = FALSE
    )
    if (isSubset(geneIDType, c("HGNC", "OMIM"))) {
        x <- as.integer(x)
    }
    df[["geneId"]] <- x
    df[["dbXrefs"]] <- NULL
    keep <- !duplicated(df[["geneId"]])
    df <- df[keep, , drop = FALSE]
    assert(hasNoDuplicates(df[["geneId"]]))
    df <- df[order(df[["geneId"]]), , drop = FALSE]
    rownames(df) <- df[["geneId"]]
    metadata(df) <- list(
        "date" = Sys.Date(),
        "geneIDType" = geneIDType,
        "organism" = organism,
        "taxonomicGroup" = taxonomicGroup
    )
    df
}
