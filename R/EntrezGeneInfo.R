#' Import NCBI Entrez gene identifier information
#'
#' @export
#' @note Updated 2021-02-12.
#'
#' @inheritParams AcidRoxygen::params
#' @inheritParams params
#'
#' @return `EntrezGeneInfo`.
#'
#' @seealso
#' - ftp://ftp.ncbi.nih.gov/gene/DATA/GENE_INFO/
#' - [geneSynonyms()].
#'
#' @examples
#' object <- EntrezGeneInfo(
#'     organism = "Homo sapiens",
#'     taxonomicGroup = "Mammalia"
#' )
#' print(object)
EntrezGeneInfo <-  # nolint
    function(
        organism,
        taxonomicGroup = NULL
    ) {
        assert(
            hasInternet(),
            isOrganism(organism),
            isString(taxonomicGroup, nullOK = TRUE)
        )
        baseURL <- pasteURL(
            "ftp.ncbi.nih.gov",
            "gene", "DATA", "GENE_INFO",
            protocol = "ftp"
        )
        if (is.null(taxonomicGroup)) {
            taxonomicGroup <- .matchNcbiTaxonomicGroup(
                organism = organism,
                mode = "geneInfo"
            )
        }
        url <- pasteURL(
            baseURL,
            taxonomicGroup,
            paste0(
                gsub(" ", "_", organism),
                ".gene_info.gz"
            )
        )
        df <- import(file = .cacheIt(url), format = "tsv", colnames = TRUE)
        df <- as(df, "DataFrame")
        colnames(df) <- camelCase(colnames(df), strict = TRUE)
        ## Ensure that any "-" values get sanitized to "NA".
        assert(any(is.na(df[["locusTag"]])))
        splitToList <- function(x) {
            x <- strsplit(x = x, split = "|", fixed = TRUE)
            x <- CharacterList(x)
            x <- sort(unique(x))
            x
        }
        df[["dbXrefs"]] <- splitToList(df[["dbXrefs"]])
        df[["synonyms"]] <- splitToList(df[["synonyms"]])
        assert(
            isSubset(
                x = c("geneId", "symbol"),
                y = colnames(df)
            ),
            hasNoDuplicates(df[["geneId"]])
        )
        colnames(df)[colnames(df) == "symbol"] <- "geneName"
        df <- removeNA(df)
        df <- df[, sort(colnames(df))]
        rownames(df) <- df[["geneId"]]
        df <- encode(df)
        metadata(df) <- list(
            "date" = Sys.Date(),
            "organism" = organism,
            "taxonomicGroup" = "taxonomicGroup"
        )
        new(Class = "EntrezGeneInfo", df)
    }
