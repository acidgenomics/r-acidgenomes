#' Import gene ontology (GO) terms per gene name
#'
#' @export
#' @note Updated 2023-12-15.
#'
#' @inheritParams AcidRoxygen::params
#'
#' @param geneNames `character` or `NULL`.
#' Return specific genes, or when `NULL`, return all genes.
#'
#' @param format `character(1)`.
#' Return GO terms in `"long"`, `"split"`, or `"nested"` formats.
#' Processing GO terms into nested format is CPU intensive, and is recommended
#' to set `geneNames` in this case.
#'
#' @return Varies, depending on `format` argument.
#'
#' - `"long"`: `DFrame` with `"geneName"`, `"goCategory"`, `"goId"`, and
#'   `"goName"` columns.
#' - `"split"`: `SplitDFrameList`: Split by `geneName` column. Contains
#'   `"geneName"`, `"goCategory"`, `"goId"`, and `"goName"` columns.
#' - `"nested"`: `DFrame` with `"geneName"`, `"goBp"`, `"goCc"`, and `"goMf"`
#'   columns. Does not contain any duplicates in `"geneName"` column.
#'
#' @seealso
#' - http://current.geneontology.org/products/pages/downloads.html
#' - https://www.ebi.ac.uk/GOA/
#' - https://www.ncbi.nlm.nih.gov/gene/
#'
#' @examples
#' object <- goTermsPerGeneName(organism = "Homo sapiens", format = "long")
#' print(object)
goTermsPerGeneName <-
    function(
        organism,
        geneNames = NULL,
        format = c("long", "split", "nested")) {
        assert(
            isCharacter(geneNames, nullOk = TRUE),
            isOrganism(organism)
        )
        organism <- match.arg(
            arg = organism,
            choices = c("Homo sapiens", "Mus musculus")
        )
        format <- match.arg(format)
        gafFile <- switch(
            EXPR = organism,
            "Homo sapiens" = "goa_human.gaf.gz",
            "Mus musculus" = "mgi.gaf.gz"
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
        gaf <- import(con = .cacheIt(url), format = "gaf")
        df <- as.data.frame(gaf)
        df <- as(df, "DFrame")
        colnames(df) <- camelCase(colnames(df))
        df <- df[, c("elements", "sets", "aspect")]
        colnames(df) <- c("geneName", "goId", "goCategory")
        df <- df[complete.cases(df), , drop = FALSE]
        if (!is.null(geneNames)) {
            i <- df[["geneName"]] %in% geneNames
            assert(any(i), msg = "Failed to match against any gene names.")
            df <- df[i, ]
        }
        df <- unique(df)
        df <- df[, c("geneName", "goCategory", "goId")]
        df <- sort(df)
        df <- leftJoin(df, goMap, by = "goId")
        if (isSubset(format, c("nested", "split"))) {
            alert("Splitting GO terms by gene.")
            spl <- split(x = df, f = df[["geneName"]])
        }
        if (identical(format, "nested")) {
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
        }
        metadata(df) <- list(
            "date" = Sys.Date(),
            "format" = format,
            "geneNames" = geneNames,
            "organism" = organism,
            "url" = url
        )
        df
    }
