## FIXME Can we add gene ontology information?
## http://current.geneontology.org/products/pages/downloads.html
## http://geneontology.org/gene-associations/goa_human_rna.gaf.gz



#' Import NCBI (Entrez) gene identifier information
#'
#' @export
#' @note Updated 2023-12-13.
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
            go <- .goTermsPerGeneName(
                organism = organism,
                geneNames = unique(df[["geneName"]])
            )
            df <- leftJoin(df, go, by = "geneName")
        }
        df <- encode(df)
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



#' Import nested GO terms per gene name
#'
#' @note Updated 2023-12-13.
#' @noRd
#'
#' @seealso
#' - http://current.geneontology.org/products/pages/downloads.html
#' - https://www.ebi.ac.uk/GOA/
#' - https://www.ncbi.nlm.nih.gov/gene/
.goTermsPerGeneName <- function(organism, geneNames) {
    assert(
        isOrganism(organism),
        isCharacter(geneNames),
        hasNoDuplicates(geneNames)
    )
    gafFile <- switch(
        EXPR = organism,
        "Homo sapiens" = "goa_human.gaf.gz",
        "Mus musculus" = "mgi.gaf.gz",
        abort(sprintf(
            "GO term lookup for {.var %s} not currently supported.",
            organism
        ))
    )
    goMap <- mapGoTerms()
    assert(identical(c("id", "name"), colnames(goMap)))
    colnames(goMap) <- c("goId", "goName")
    url <- pasteUrl(
        "geneontology.org",
        "gene-associations",
        gafFile,
        protocol = "https"
    )
    gaf <- import(
        con = .cacheIt(url),
        format = "gaf"
    )
    df <- as.data.frame(gaf)
    df <- as(df, "DFrame")
    colnames(df) <- camelCase(colnames(df))
    df <- df[, c("elements", "sets", "aspect")]
    colnames(df) <- c("geneName", "goId", "goCategory")
    df <- df[complete.cases(df), ]
    i <- df[["geneName"]] %in% geneNames
    assert(any(i), msg = "Failed to match against any gene names.")
    df <- df[i, ]
    df <- unique(df)
    df <- df[, c("geneName", "goCategory", "goId")]
    df <- sort(df)
    df <- leftJoin(df, goMap, by = "goId")
    spl <- split(x = df, f = df[["geneName"]])
    alert("Nesting GO terms per gene.")
    lst <- mclapply(
        X = spl,
        FUN = function(x) {
            idx <- list(
                "bp" = which(x[["goCategory"]] == "BP"),
                "cc" = which(x[["goCategory"]] == "CC"),
                "mf" = which(x[["goCategory"]] == "MF")
            )
            out <- DataFrame(
                "geneName" = x[["geneName"]][[1L]],
                "goBp" = I(list(paste(
                    x[["goId"]][idx[["bp"]]],
                    x[["goName"]][idx[["bp"]]]
                ))),
                "goCc" = I(list(paste(
                    x[["goId"]][idx[["cc"]]],
                    x[["goName"]][idx[["cc"]]]
                ))),
                "goMf" = I(list(paste(
                    x[["goId"]][idx[["mf"]]],
                    x[["goName"]][idx[["mf"]]]
                )))
            )
            out
        }
    )
    dfl <- DataFrameList(lst)
    df <- unlist(dfl)
    assert(is(df, "DFrame"))
    df[["goBp"]] <- CharacterList(df[["goBp"]])
    df[["goCc"]] <- CharacterList(df[["goCc"]])
    df[["goMf"]] <- CharacterList(df[["goMf"]])
    metadata(df) <- list(
        "date" = Sys.Date(),
        "geneNames" = geneNames,
        "organism" = organism,
        "url" = url
    )
    df
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
