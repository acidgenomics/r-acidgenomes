#' Import NCBI (Entrez) gene identifier information
#'
#' @export
#' @note Updated 2023-12-12.
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
            isString(taxonomicGroup, nullOk = TRUE),
            isFlag(cache)
        )
        baseURL <- pasteUrl(
            "ftp.ncbi.nih.gov", "gene", "DATA", "GENE_INFO",
            protocol = "https"
        )
        if (is.null(taxonomicGroup)) {
            taxonomicGroup <- .matchNcbiTaxonomicGroup(
                organism = organism,
                mode = "geneInfo"
            )
        }
        url <- pasteUrl(
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
        df <- import(
            con = ifelse(
                test = cache,
                yes = .cacheIt(url),
                no = url
            ),
            format = "tsv",
            colnames = TRUE,
            naStrings = "-"
        )
        df <- as(df, "DFrame")
        colnames(df) <- camelCase(colnames(df), strict = TRUE)
        assert(
            isSubset(
                x = c(
                    "geneId",
                    "locusTag",
                    "modificationDate",
                    "symbol",
                    "synonyms",
                    "xTaxId"
                ),
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
        df <- removeNa(df)
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
        df[["modificationDate"]] <- sub(
            pattern = "^([0-9]{4})([0-9]{2})([0-9]{2})$",
            replacement = "\\1-\\2-\\3",
            x = df[["modificationDate"]]
        )
        df[["modificationDate"]] <- as.Date(df[["modificationDate"]])
        if (isInstalled("GeneSummary")) {
            gs <- .refseqGeneSummary(organism)
            df <- leftJoin(df, gs, by = "geneId")
            metadata(df)[["refseqGeneSummary"]] <- TRUE
        }
        df <- encode(df)
        metadata(df)[["date"]] <- Sys.Date()
        metadata(df)[["organism"]] <- organism
        metadata(df)[["packageVersion"]] <- .pkgVersion
        metadata(df)[["taxonomicGroup"]] <- taxonomicGroup
        metadata(df)[["url"]] <- url
        new(Class = "NcbiGeneInfo", df)
    }



#' Get RefSeq gene summary
#'
#' @note Updated 2023-12-12.
#' @noRd
#'
#' @seealso
#' - https://www.ncbi.nlm.nih.gov/refseq/about/
#' - https://ftp.ncbi.nih.gov/refseq/release/complete/
.refseqGeneSummary <- function(organism) {
    assert(requireNamespaces("GeneSummary"))
    taxId <- .mapOrganismToNcbiTaxId(organism)
    df <- GeneSummary::loadGeneSummary(organism = taxId)
    assert(is.data.frame(df))
    df <- as(df, "DFrame")
    colnames(df) <- camelCase(colnames(df))
    cols <- c(
        ## > "refSeqAccession",
        ## > "organism",
        ## > "taxonId",
        "geneId",
        "reviewStatus",
        "geneSummary"
    )
    assert(isSubset(cols, colnames(df)))
    df <- df[, cols, drop = FALSE]
    df <- df[complete.cases(df), , drop = FALSE]
    df <- unique(df)
    df[["reviewStatus"]] <- factor(
        x = df[["reviewStatus"]],
        levels = c(
            "REVIEWED REFSEQ",
            "VALIDATED REFSEQ",
            "PROVISIONAL REFSEQ",
            "INFERRED REFSEQ",
            "PREDICTED REFSEQ"
        )
    )
    df <- df[order(df), , drop = FALSE]
    df <- df[!duplicated(df[["geneId"]]), , drop = FALSE]
    assert(hasNoDuplicates(df[["geneId"]]))
    df[["reviewStatus"]] <- NULL
    df
}
