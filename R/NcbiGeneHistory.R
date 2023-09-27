#' NCBI gene history
#'
#' @export
#' @note Updated 2023-09-26.
#'
#' @inheritParams AcidRoxygen::params
#'
#' @return `NcbiGeneHistory`.
#'
#' @examples
#' ## Homo sapiens.
#' object <- NcbiGeneHistory(organism = "Homo sapiens")
#' print(object)
NcbiGeneHistory <- # nolint
    function(organism) {
        assert(isOrganism(organism))
        taxId <- .mapOrganismToNcbiTaxId(organism)
        assert(isInt(taxId))
        url <- pasteURL(
            "ftp.ncbi.nih.gov", "gene", "DATA", "gene_history.gz",
            protocol = "https"
        )
        ## readr is much faster than base engine at parsing this file.
        df <- import(
            con = .cacheIt(url),
            format = "tsv",
            engine = ifelse(
                test = isInstalled("readr"),
                yes = "readr",
                no = "base"
            ),
            naStrings = "-"
        )
        df <- as(df, "DFrame")
        colnames(df) <- camelCase(colnames(df))
        assert(identical(
            x = colnames(df),
            y = c(
                "xTaxId",
                "geneId",
                "discontinuedGeneId",
                "discontinuedSymbol",
                "discontinueDate"
            ),
            hasNoDuplicates(df[["discontinuedGeneId"]])
        ))
        keep <- df[["xTaxId"]] == taxId
        df[["xTaxId"]] <- NULL
        df <- df[keep, , drop = FALSE]
        df[["discontinueDate"]] <- sub(
            pattern = "^([0-9]{4})([0-9]{2})([0-9]{2})$",
            replacement = "\\1-\\2-\\3",
            x = df[["discontinueDate"]]
        )
        df[["discontinueDate"]] <- as.Date(df[["discontinueDate"]])
        df[["discontinuedGeneId"]] <- as.integer(df[["discontinuedGeneId"]])
        df[["geneId"]] <- as.integer(df[["geneId"]])
        j <- c(
            "discontinuedGeneId",
            "discontinuedSymbol",
            "discontinueDate",
            "geneId"
        )
        df <- df[, j, drop = FALSE]
        i <- order(df)
        df <- df[i, , drop = FALSE]
        df <- encode(df)
        rownames(df) <- df[["discontinuedGeneId"]]
        metadata(df) <- list(
            "date" = Sys.Date(),
            "organism" = organism,
            "packageVersion" = .pkgVersion,
            "taxonomyId" = taxId,
            "url" = url
        )
        new(df, Class = "NcbiGeneHistory")
    }
