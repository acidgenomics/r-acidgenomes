#' Import NCBI (Entrez) gene identifier information
#'
#' @export
#' @note Updated 2023-12-19.
#'
#' @inheritParams AcidRoxygen::params
#'
#' @param taxonomicGroup `character(1)`.
#' NCBI FTP server taxonomic group subdirectory path (e.g. "Mammalia").
#' Defining this manually avoids having to query the FTP server.
#'
#' @param refseqGeneSummary `logical(1)`.
#' Include RefSeq gene summary in `"refseqGeneSummary"` column.
#' Requires Bioconductor GeneSummary package to be installed.
#'
#' @param goTerms `logical(1)`.
#' Return nested gene ontology (GO) terms in `"goBp"` (BP: biological process),
#' `"goCc"` (CC: cellular component), and `"goMf"` (MF: molecular function)
#' columns. This is computationally intensive, and not supported for all
#' genomes. Intended primarily for *Homo sapiens*.
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
             refseqGeneSummary = FALSE,
             goTerms = FALSE) {
        assert(
            hasInternet(),
            isOrganism(organism),
            isString(taxonomicGroup, nullOk = TRUE),
            isFlag(refseqGeneSummary),
            isFlag(goTerms)
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
            con = .cacheIt(url),
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
        if (isTRUE(refseqGeneSummary)) {
            gs <- .refseqGeneSummary(organism)
            df <- leftJoin(df, gs, by = "geneId")
        }
        if (isTRUE(goTerms)) {
            go <- goTermsPerGeneName(
                organism = organism,
                geneNames = unique(df[["geneName"]]),
                format = "nested"
            )
            df <- leftJoin(df, go, by = "geneName")
        }
        ## Disabled Rle encoding in 0.7.3 update.
        ## > df <- encode(df)
        df <- df[, sort(colnames(df))]
        metadata(df) <- list(
            "date" = Sys.Date(),
            "organism" = organism,
            "packageVersion" = .pkgVersion,
            "refseqGeneSummary" = refseqGeneSummary,
            "taxonomicGroup" = taxonomicGroup,
            "url" = url
        )
        new(Class = "NcbiGeneInfo", df)
    }



#' Import RefSeq gene summary
#'
#' @note Updated 2023-12-12.
#' @noRd
#'
#' @seealso
#' - https://www.ncbi.nlm.nih.gov/refseq/about/
#' - https://ftp.ncbi.nih.gov/refseq/release/complete/
.refseqGeneSummary <- function(organism) {
    assert(
        requireNamespaces("GeneSummary"),
        isOrganism(organism)
    )
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
    colnames(df)[colnames(df) == "geneSummary"] <- "refseqGeneSummary"
    df[["reviewStatus"]] <- NULL
    df
}
