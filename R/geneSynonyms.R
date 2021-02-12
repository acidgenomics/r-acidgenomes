## FIXME NEED TO SPECIFY AN ENSEMBL MODE HERE.
## FIXME CONSIDER USING THIS METADATA FOR DEPMAPANALYSIS...
## FIXME This can contain some Ensembl duplicates (e.g. ENSG00000004866).



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
#' ## > print(object)
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
    df <- df[, cols]
    keep <- !all(is.na(df[["geneSynonyms"]]))
    df <- df[keep, , drop = FALSE]
    if (identical(geneIDType, "Entrez")) {
        df[["dbXrefs"]] <- NULL
        return(df)
    }
    pattern <- paste0(
        "^",
        switch(
            EXPR = geneIDType,
            "OMIM" = "MIM",
            geneIDType
        ),
        ":(.+)$"
    )
    keep <- any(grepl(pattern = pattern, x = df[["dbXrefs"]]))
    df <- df[keep, , drop = FALSE]



    ## FIXME NEED TO REWORK LIST APPROACH HERE...LAPPLY?
    ## FIXME NEED TO EXTRACT THE ELEMENTS THAT MATCH AND MAKE IT FLAG....
    xxx <- gsub(
        pattern = pattern,
        replacement = "\\1",
        x = df[["dbXrefs"]][lgl]
    )
    xxx <- unlist(xxx, recursive = FALSE, use.names = FALSE)
    assert(identical(length(xxx), nrow(df)))

    length(xxx)
    ## FIXME DONT DO THIS, AS IT WILL DROP IDENTIFIERS.
    yyy <- unlist(xxx, recursive = FALSE, use.names = TRUE)
    length(yyy)









    df <- df[, c("symbol", "synonyms", "dbXrefs")]
    colnames(df)[colnames(df) == "symbol"] <- "geneName"
    colnames(df)[colnames(df) == "synonyms"] <- "geneSynonyms"
    keep <- df[["geneSynonyms"]] != "-"
    df <- df[keep, , drop = FALSE]
    keep <- df[["dbXrefs"]] != "-"
    df <- df[keep, , drop = FALSE]



    ## FIXME NEED TO IMPROVE SUPPORT FOR ID TYPE RETURN HERE.
    ## FIXME ENSEMBL MATCH MIGHT NOT WORK FOR OTHER GENOMES....


    ## Sanitize the identifiers.
    pattern <- switch(
        EXPR = organism,
        "Drosophila melanogaster" = {
            "\\bFBgn[0-9]{7}\\b"
        },
        "\\bENS[A-Z]+[0-9]{11}\\b"
    )
    df[["geneId"]] <- str_extract(
        string = df[["dbXrefs"]],
        pattern = pattern
    )
    keep <- !is.na(df[["geneId"]])
    df <- df[keep, , drop = FALSE]
    cols <- c("geneId", "geneSynonyms")
    df <- df[, cols, drop = FALSE]
    df <- df[order(df[[1L]]), , drop = FALSE]
    x <- split(df, f = df[[1L]])
    if (identical(return, "SplitDataFrameList")) {
        return(x)
    }
    alert("Preparing unique synonyms per gene.")


    xx <- lapply(
        X = x[, 2L],
        FUN = function(x) {
            geneId <- x[["geneId"]][[1L]]
            geneSynonyms <- strsplit(
                x = sort(
                    x = x[["geneSynonyms"]],
                    decreasing = FALSE,
                    na.last = TRUE
                ),
                split = ", "
            )[[1L]]
            geneSynonyms <- unlist(geneSynonyms)
            geneSynonyms <- c(geneSynonyms, x[["geneName"]])
            geneSynonyms <- sort(unique(geneSynonyms))
            geneSynonyms <- toString(geneSynonyms)
            DataFrame("geneId" = geneId, "geneSynonyms" = geneSynonyms)
        }
    )



    df <- do.call(what = rbind, args = list)
    assert(identical(names(split), df[["geneId"]]))
    df <- df[complete.cases(df), ]
    rownames(df) <- df[["geneId"]]
    out <- df
    out
}
