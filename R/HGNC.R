#' Import Human Genome Organization (HUGO) Gene Nomenclature Committee (HGNC)
#' metadata
#'
#' @export
#' @note Updated 2023-09-27.
#'
#' @return `HGNC`.
#'
#' @seealso
#' - https://www.genenames.org/
#' - https://www.genenames.org/download/statistics-and-files/
#'
#' @examples
#' object <- HGNC()
#' print(object)
HGNC <- # nolint
    function() {
        alert("Importing HGNC complete set.")
        url <- pasteURL(
            "ftp.ebi.ac.uk",
            "pub",
            "databases",
            "genenames",
            "new",
            "tsv",
            "hgnc_complete_set.txt",
            protocol = "https"
        )
        file <- .cacheIt(url)
        lines <- import(file, format = "lines")
        ## FIXME Rework this as `fillLines` in pipette and then use here.
        spl <- strsplit(x = lines, split = "\t", fixed = TRUE)
        fixIdx <- which(lengths(spl) != length(spl[[1L]]))
        spl[fixIdx] <- lapply(
            X = spl[fixIdx],
            FUN = function(x) {
                append(x = x, values = NA)
            }
        )
        assert(all(lengths(spl) == length(spl[[1L]])))
        lines <- vapply(
            X = spl,
            FUN = paste0,
            collapse = "\t",
            FUN.VALUE = character(1L)
        )
        con <- textConnection(lines)
        df <- import(con = con, format = "tsv")
        close(con)
        df <- as(df, "DFrame")
        colnames(df) <- camelCase(colnames(df), strict = TRUE)
        df[["agr"]] <- NULL
        df[["gencc"]] <- NULL
        assert(
            isSubset(
                x = c(
                    "dateApprovedReserved",
                    "dateModified",
                    "dateNameChanged",
                    "dateSymbolChanged",
                    "entrezId",
                    "hgncId"
                ),
                y = colnames(df)
            ),
            hasNoDuplicates(df[["hgncId"]])
        )
        df[["hgncId"]] <- as.integer(sub(
            pattern = "^HGNC\\:",
            replacement = "",
            x = df[["hgncId"]]
        ))
        df <- df[order(df[["hgncId"]]), , drop = FALSE]
        rownames(df) <- df[["hgncId"]]
        isNested <- bapply(
            X = df,
            FUN = function(x) {
                any(grepl(pattern = "|", x = x, fixed = TRUE))
            }
        )
        if (any(isNested)) {
            vars <- names(isNested)[isNested]
            df <- mutateAt(
                object = df,
                vars = vars,
                fun = .splitToCharacterList,
                split = "|"
            )
        }
        colnames(df)[colnames(df) == "entrezId"] <- "ncbiGeneId"
        df[["dateApprovedReserved"]] <- as.Date(df[["dateApprovedReserved"]])
        df[["dateModified"]] <- as.Date(df[["dateModified"]])
        df[["dateNameChanged"]] <- as.Date(df[["dateNameChanged"]])
        df[["dateSymbolChanged"]] <- as.Date(df[["dateSymbolChanged"]])
        df[["ncbiGeneId"]] <- as.integer(df[["ncbiGeneId"]])
        metadata(df) <- list(
            "date" = Sys.Date(),
            "packageVersion" = .pkgVersion,
            "url" = url
        )
        new("HGNC", df)
    }



## Updated 2021-03-19.
.splitToCharacterList <- function(x, split = "|") {
    if (all(is.na(x))) {
        return(x)
    }
    x <- strsplit(x, split = split, fixed = TRUE)
    x <- CharacterList(x)
    ## > x <- sort(unique(x))
    x
}
