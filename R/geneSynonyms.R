## FIXME NEED TO SPECIFY AN ENSEMBL MODE HERE.
## FIXME CONSIDER USING THIS METADATA FOR DEPMAPANALYSIS...



## NOTE This can contain some Ensembl duplicates (e.g. ENSG00000004866).

## FIXME CONSIDER ADDING SUPPORT FOR THESE ORGANISMS:
## Bos_taurus.gene_info.gz	1.6 MB	2/9/21, 10:44:00 PM
## Canis_familiaris.gene_info.gz	1.3 MB	2/9/21, 10:44:00 PM
## Homo_sapiens.gene_info.gz	2.8 MB	2/9/21, 10:44:00 PM
## Mus_musculus.gene_info.gz	3.0 MB	2/9/21, 10:44:00 PM
## Pan_troglodytes.gene_info.gz	1.3 MB	2/9/21, 10:44:00 PM
## Rattus_norvegicus.gene_info.gz	2.1 MB	2/9/21, 10:44:00 PM
## Sus_scrofa.gene_info.gz	1.2 MB	2/9/21, 10:44:00 PM



#' Gene synonyms
#'
#' Look up gene synonyms from NCBI.
#'
#' @note Updated 2021-02-10.
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
#' ## > object <- geneSynonyms(organism = "Homo sapiens")
#' ## > print(object)
geneSynonyms <- function(
    organism = c(
        "Homo sapiens",
        "Mus musculus",
        "Drosophila melanogaster"
    ),
    idType = c("Entrez", "Ensembl", "HGNC", "OMIM")
) {
    organism <- match.arg(organism)
    idType <- match.arg(idType)
    assert(
        hasInternet(),
        isOrganism(organism)
    )
    alert(sprintf(
        "Importing {.var %s} synonyms from NCBI for %s gene identifiers.",
        organism, idType
    ))
    genome <- c(
        "kingdom" = switch(
            EXPR = species,
            "Drosophila melanogaster" = "Invertebrates",
            "Mammalia"
        ),
        ## NCBI uses underscore for species name.
        "species" = gsub(" ", "_", organism)
    )
    url <- pasteURL(
        "ftp.ncbi.nih.gov",
        "gene",
        "DATA",
        "GENE_INFO",
        genome[["kingdom"]],
        paste0(genome[["species"]], ".gene_info.gz"),
        protocol = "ftp"
    )
    df <- import(file = .cacheIt(url), format = "tsv", colnames = TRUE)
    df <- as(df, "DataFrame")
    colnames(df) <- camelCase(colnames(df), strict = TRUE)
    df <- df[, c("geneId", "synonyms", "dbXrefs")]
    keep <- df[["synonyms"]] != "-"
    df <- df[keep, ]
    ## Ensure synonyms include current gene symbol.
    df[["synonyms"]] <- paste(df[["synonyms"]], df[["symbol"]], sep = "|")
    df[["symbol"]] <- NULL
    splitToList <- function(x) {
        x <- strsplit(x = x, split = "|", fixed = TRUE)
        x <- CharacterList(x)
        x <- sort(unique(x))
        x
    }
    df[["synonyms"]] <- splitToList(df[["synonyms"]])
    colnames(df)[colnames(df) == "synonyms"] <- "geneSynonyms"
    if (identical(idType, "Entrez")) {
        df[["dbXrefs"]] <- NULL
        return(df)
    }
    df[["dbXrefs"]] <- splitToList(df[["dbXrefs"]])
    pattern <- paste0(
        "^",
        switch(
            EXPR = idType,
            "OMIM" = "MIM",
            idType
        ),
        ":(.+)$"
    )
    keep <- any(grepl(pattern = pattern, x = df[["dbXrefs"]]))
    df <- df[keep, , drop = FALSE]
    lgl <- grepl(pattern = pattern, x = df[["dbXrefs"]])
    assert(is(lgl, "LogicalList"))

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
