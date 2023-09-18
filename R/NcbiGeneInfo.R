## FIXME Ensure modification date is classed as a date.



#' Import NCBI (Entrez) gene identifier information
#'
#' @export
#' @note Updated 2023-09-15.
#'
#' @inheritParams AcidRoxygen::params
#'
#' @param cache `logical(1)`.
#' Cache the gene info file from NCBI FTP server using BiocFileCache.
#'
#' @param taxonomicGroup `character(1)`.
#' NCBI FTP server taxonomic group subdirectory path (e.g. "Mammalia").
#' Defining this manually avoids having to query the FTP server.
#'
#' @return `NcbiGeneInfo`.
#'
#' @seealso
#' - https://ftp.ncbi.nih.gov/gene/DATA/GENE_INFO/
#'
#' @examples
#' object <- NcbiGeneInfo(
#'     organism = "Homo sapiens",
#'     taxonomicGroup = "Mammalia"
#' )
#' print(object)
NcbiGeneInfo <- # nolint
    function(organism,
             taxonomicGroup = NULL,
             cache = TRUE) {
        assert(
            hasInternet(),
            isOrganism(organism),
            isString(taxonomicGroup, nullOK = TRUE),
            isFlag(cache)
        )
        baseURL <- pasteURL(
            "ftp.ncbi.nih.gov", "gene", "DATA", "GENE_INFO",
            protocol = "https"
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
        df <- import(
            con = ifelse(
                test = cache,
                yes = .cacheIt(url),
                no = url
            ),
            format = "tsv",
            colnames = TRUE,
            engine = "readr"
        )
        df <- as(df, "DFrame")
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
        colnames(df)[colnames(df) == "xTaxId"] <- "taxonomyId"
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
        ## FIXME Include package name and package version.
        metadata(df) <- list(
            "date" = Sys.Date(),
            "organism" = organism,
            "taxonomicGroup" = "taxonomicGroup",
            "url" = url
        )
        new(Class = "NcbiGeneInfo", df)
    }
