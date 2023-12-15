## FIXME Can't process Mus musculus currently.
## mgi.gaf.gz can't import into a TidySet.



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
#' - https://geneontology.org/docs/go-annotation-file-gaf-format-2.2/
#' - https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2go.gz
#' - Python goatools package.
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
        format <- match.arg(format)
        gafFile <- switch(
            EXPR = organism,
            "Arabidopsis thaliana" = "tair.gaf.gz",
            "Bos taurus" = "goa_cow.gaf.gz",
            "Caenorhabditis elegans" = "wb.gaf.gz",
            "Canis lupus familiaris" = "goa_dog.gaf.gz",
            "Danio rerio" = "zfin.gaf",
            "Drosophila melanogaster" = "fb.gaf.gz",
            "Gallus gallus" = "goa_chicken.gaf.gz",
            "Homo sapiens" = "goa_human.gaf.gz",
            "Mus musculus" = "mgi.gaf.gz",
            "Rattus norvegicus" = "rgd.gaf.gz",
            "Saccharomyces cerevisiae" = "sgd.gaf.gz",
            "Sus scrofa" = "goa_pig.gaf.gz",
            abort(sprintf("Unsupported organism: {.var %s}.", organism))
        )
        url <- pasteUrl(
            "geneontology.org",
            "gene-associations",
            gafFile,
            protocol = "https"
        )
        df <- import(con = .cacheIt(url), format = "gaf")
        assert(is.data.frame(df))
        df <- as(df, "DFrame")
        cols <- c("dbObjectSymbol", "goId", "aspect")
        assert(isSubset(cols, colnames(df)))
        df <- df[, cols]
        colnames(df) <- c("geneName", "goId", "goCategory")
        df <- df[complete.cases(df), ]
        assert(identical(
            x = sort(unique(df[["goCategory"]])),
            y = c("C", "F", "P")
        ))
        df[["goCategory"]] <- sub(
            pattern = "C",
            replacement = "CC",
            x = df[["goCategory"]],
            fixed = TRUE
        )
        df[["goCategory"]] <- sub(
            pattern = "F",
            replacement = "MF",
            x = df[["goCategory"]],
            fixed = TRUE
        )
        df[["goCategory"]] <- sub(
            pattern = "P",
            replacement = "BP",
            x = df[["goCategory"]],
            fixed = TRUE
        )
        if (!is.null(geneNames)) {
            i <- df[["geneName"]] %in% geneNames
            assert(any(i), msg = "Failed to match against any gene names.")
            df <- df[i, ]
        }
        df <- unique(df)
        df <- df[, c("geneName", "goCategory", "goId")]
        df <- sort(df)
        goMap <- mapGoTerms()
        assert(identical(c("id", "name"), colnames(goMap)))
        colnames(goMap) <- c("goId", "goName")
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
