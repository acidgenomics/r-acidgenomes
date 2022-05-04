#' Gene synonyms
#'
#' Look up gene synonyms from NCBI.
#'
#' @note Updated 2022-05-04.
#' @export
#'
#' @section *Caenorhabditis elegans*:
#'
#' Synonym support for *Caenorhabditis elegans* is poor on NCBI. Use the
#' [WormBase](https://r.acidgenomics.com/packages/wormbase/) package instead.
#'
#'
#' @inheritParams AcidRoxygen::params
#' @inheritParams EntrezGeneInfo
#' @inheritParams params
#' @param geneIdType `character(1)`.
#' Type of gene identifier to return in the `geneId` column.
#'
#' @return `DataFrame` containing `geneId` and `geneSynonyms` columns.
#'
#' @examples
#' object <- geneSynonyms(
#'     organism = "Homo sapiens",
#'     taxonomicGroup = "Mammalia",
#'     geneIdType = "Ensembl"
#' )
#' print(object)
geneSynonyms <-
    function(organism,
             taxonomicGroup = NULL,
             geneIdType = c("Entrez", "Ensembl", "HGNC", "OMIM")) {
        geneIdType <- match.arg(geneIdType)
        df <- EntrezGeneInfo(
            organism = organism,
            taxonomicGroup = taxonomicGroup
        )
        meta <- metadata(df)
        meta <- append(x = meta, values = list("geneIdType" = geneIdType))
        cols <- c("geneId", "geneSynonyms", "dbXrefs")
        assert(
            isSubset(cols, colnames(df)),
            is(df[["geneSynonyms"]], "CharacterList")
        )
        df <- as(df, "DataFrame")
        df <- df[, cols, drop = FALSE]
        keep <- !all(is.na(df[["geneSynonyms"]]))
        df <- df[keep, , drop = FALSE]
        if (identical(geneIdType, "Entrez")) {
            df[["geneId"]] <- decode(df[["geneId"]])
            assert(is.integer(df[["geneId"]]))
            df[["dbXrefs"]] <- NULL
            return(df)
        }
        pattern <- paste0(
            "^",
            switch(
                EXPR = geneIdType,
                "HGNC" = "HGNC:HGNC",
                "OMIM" = "MIM",
                geneIdType
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
        if (isSubset(geneIdType, c("HGNC", "OMIM"))) {
            x <- as.integer(x)
        }
        df[["geneId"]] <- x
        df[["dbXrefs"]] <- NULL
        keep <- !duplicated(df[["geneId"]])
        df <- df[keep, , drop = FALSE]
        assert(hasNoDuplicates(df[["geneId"]]))
        df <- df[order(df[["geneId"]]), , drop = FALSE]
        rownames(df) <- df[["geneId"]]
        metadata(df) <- meta
        df
    }
