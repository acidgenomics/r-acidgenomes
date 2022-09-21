## FIXME Add `cache = TRUE` to function here.
## FIXME Need to rename `xTaxId` to `taxId`.



#' Import NCBI Entrez gene identifier information
#'
#' @export
#' @note Updated 2022-09-21.
#'
#' @inheritParams AcidRoxygen::params
#' @param taxonomicGroup `character(1)`.
#' NCBI FTP server taxonomic group subdirectory path (e.g. "Mammalia").
#' Defining this manually avoids having to query the FTP server.
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
EntrezGeneInfo <- # nolint
    function(organism,
             taxonomicGroup = NULL) {
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
        alert(sprintf(
            "Downloading {.emph %s} gene info from NCBI at {.url %s}.",
            organism, url
        ))
        ## Input TSV is malformed, as of 2022-05-04. Readr handles this more
        ## gracefully than the base import engine.
        ## Error: line 47 did not have 16 elements
        df <- import(
            con = .cacheIt(url),
            format = "tsv",
            colnames = TRUE,
            engine = "readr"
        )
        df <- as(df, "DataFrame")
        colnames(df) <- camelCase(colnames(df), strict = TRUE)
        assert(
            isSubset(
                x = c("geneId", "locusTag", "symbol", "synonyms"),
                y = colnames(df)
            ),
            hasNoDuplicates(df[["geneId"]]),
            ## Ensure that any "-" values get sanitized to "NA".
            anyNA(df[["locusTag"]])
        )
        df[["fullNameFromNomenclatureAuthority"]] <- NULL
        df[["symbolFromNomenclatureAuthority"]] <- NULL
        df[["geneId"]] <- as.integer(df[["geneId"]])
        colnames(df)[colnames(df) == "symbol"] <- "geneName"
        colnames(df)[colnames(df) == "synonyms"] <- "geneSynonyms"
        df <- removeNA(df)
        df <- df[, sort(colnames(df))]
        rownames(df) <- df[["geneId"]]
        splitToList <- function(x) {
            x <- strsplit(x = x, split = "|", fixed = TRUE)
            x <- CharacterList(x)
            x <- sort(unique(x))
            x
        }
        if (isSubset("dbXrefs", colnames(df))) {
            df[["dbXrefs"]] <- splitToList(df[["dbXrefs"]])
        }
        if (isSubset("geneSynonyms", colnames(df))) {
            df[["geneSynonyms"]] <- splitToList(df[["geneSynonyms"]])
        }
        if (isSubset("otherDesignations", colnames(df))) {
            df[["otherDesignations"]] <- splitToList(df[["otherDesignations"]])
        }
        df <- encode(df)
        metadata(df) <- list(
            "date" = Sys.Date(),
            "organism" = organism,
            "taxonomicGroup" = "taxonomicGroup",
            "url" = url
        )
        new(Class = "EntrezGeneInfo", df)
    }
