#' Gene synonyms
#'
#' Look up gene synonyms from NCBI.
#'
#' @note Synonym support for *Caenorhabditis elegans* is poor on NCBI.
#' Use the [WormBase](https://r.acidgenomics.com/packages/wormbase/) package
#' instead.
#' @note Updated 2021-01-14.
#' @export
#'
#' @inheritParams AcidRoxygen::params
#'
#' @return
#' - `DataFrame`:
#'   Returns unique row for each `geneId`.
#'   Returns only `geneId` and `geneSynonyms` columns.
#'   The `geneSynonyms` column returns as character vector, with synonyms
#'   arranged alphabetically and delimited by `", "`.
#' - `SplitDataFrameList`:
#'   Split by `geneId` column.
#'   Keeps duplicate rows mapping to same `geneId` (e.g. ENSG00000004866).
#'   Returns `geneId`, `geneName`, and `geneSynonyms` columns in the split.
#'
#' @examples
#' ## CPU intensive.
#' ## > object <- geneSynonyms(organism = "Homo sapiens")
#' ## > print(object)
geneSynonyms <- function(
    organism = c(
        "Homo sapiens",
        "Mus musculus",
        "Drosophila melanogaster"
    ),
    return = c("DataFrame", "SplitDataFrameList")
) {
    assert(hasInternet())
    organism <- match.arg(organism)
    return <- match.arg(return)
    alert(sprintf("Importing {.var %s} gene synonyms from Ensembl.", organism))
    ## NCBI uses underscore for species name.
    species <- gsub(" ", "_", organism)
    if (species == "Drosophila_melanogaster") {
        ## This is covered in full local tests.
        kingdom <- "Invertebrates"  # nocov
    } else {
        kingdom <- "Mammalia"
    }
    genome <- c(kingdom = kingdom, species = species)
    url <- pasteURL(
        "ftp.ncbi.nih.gov",
        "gene",
        "DATA",
        "GENE_INFO",
        genome[["kingdom"]],
        paste0(genome[["species"]], ".gene_info.gz"),
        protocol = "ftp"
    )
    file <- .cacheIt(url)
    df <- import(file = file, format = "tsv", colnames = TRUE)
    assert(hasLength(df))
    df <- as(df, "DataFrame")
    colnames(df) <- camelCase(colnames(df), strict = TRUE)
    df <- df[, c("symbol", "synonyms", "dbXrefs")]
    colnames(df)[colnames(df) == "symbol"] <- "geneName"
    colnames(df)[colnames(df) == "synonyms"] <- "geneSynonyms"
    keep <- df[["geneSynonyms"]] != "-"
    df <- df[keep, , drop = FALSE]
    keep <- df[["dbXrefs"]] != "-"
    df <- df[keep, , drop = FALSE]
    df[["geneSynonyms"]] <- gsub(
        pattern = "\\|",
        replacement = ", ",
        x = df[["geneSynonyms"]]
    )
    ## Sanitize the identifiers.
    if (identical(organism, "Drosophila melanogaster")) {
        ## This is covered in full local tests.
        ## nocov start
        pattern <- "\\bFBgn[0-9]{7}\\b"
        ## nocov end
    } else {
        pattern <- "\\bENS[A-Z]+[0-9]{11}\\b"
    }
    df[["geneId"]] <- str_extract(
        string = df[["dbXrefs"]],
        pattern = pattern
    )
    keep <- !is.na(df[["geneId"]])
    df <- df[keep, , drop = FALSE]
    df <- df[, c("geneId", "geneName", "geneSynonyms")]
    df <- df[order(df[["geneId"]]), , drop = FALSE]
    split <- split(df, f = df[["geneId"]])
    if (identical(return, "DataFrame")) {
        alert("Preparing unique synonyms per gene.")
        list <- bplapply(
            X = split,
            FUN = function(x) {
                geneId <- x[["geneId"]][[1L]]
                geneSynonyms <- strsplit(x[["geneSynonyms"]], split = ", ")
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
    } else {
        out <- split
    }
    out
}
