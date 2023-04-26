#' Connect to AnnotationHub
#'
#' @note Updated 2023-04-26.
#' @noRd
#'
#' @details
#' On a fresh install this will print a txProgressBar to the console. We're
#' using `capture.output()` here to suppress the console output, since it's not
#' very informative and can cluster R Markdown reports.
.annotationHub <- function() {
    assert(hasInternet())
    .suppressAll({
        requireNamespaces("AnnotationHub")
        ah <- AnnotationHub::AnnotationHub(ask = FALSE)
    })
    assert(is(ah, "AnnotationHub"))
    ah
}



#' Get Ensembl / NCBI (Entrez) mappings from NCBI OrgDb via AnnotationHub
#'
#' @note Updated 2023-04-26.
#' @noRd
.getEnsembl2NcbiFromOrgDb <-
    function(keys,
             keytype,
             columns,
             organism,
             strict = TRUE) {
        .suppressAll({
            requireNamespaces(c("AnnotationDbi", "AnnotationHub"))
        })
        assert(
            isCharacter(keys),
            hasNoDuplicates(keys),
            isString(keytype),
            isCharacter(columns),
            isOrganism(organism),
            isFlag(strict)
        )
        alert(sprintf(
            "Matching identifiers using NCBI {.cls %s} via {.pkg %s} %s.",
            "OrgDb",
            "AnnotationHub",
            packageVersion("AnnotationHub")
        ))
        ah <- .annotationHub()
        ahs <- AnnotationHub::query(ah, pattern = c(organism, "NCBI", "OrgDb"))
        id <- tail(names(ahs), n = 1L)
        suppressMessages({
            orgdb <- ah[[id]]
        })
        assert(is(orgdb, "OrgDb"))
        alertInfo(sprintf(
            "{.val %s} (%s): %s.",
            id,
            mcols(ahs)[id, "title"],
            mcols(ahs)[id, "description"]
        ))
        .suppressAll({
            df <- AnnotationDbi::select(
                x = orgdb,
                keys = keys,
                keytype = keytype,
                columns = columns
            )
        })
        assert(is.data.frame(df))
        if (isTRUE(strict)) {
            df <- df[complete.cases(df), , drop = FALSE]
            if (!areSetEqual(keys, unique(df[[keytype]]))) {
                setdiff <- setdiff(keys, unique(df[[keytype]]))
                abort(sprintf(
                    "Match failure: %s.",
                    toInlineString(setdiff, n = 10L)
                ))
            }
        }
        colnames(df)[colnames(df) == "ENSEMBL"] <- "ensemblGeneId"
        colnames(df)[colnames(df) == "ENTREZID"] <- "ncbiGeneId"
        df[["ncbiGeneId"]] <- as.integer(df[["ncbiGeneId"]])
        df <- as(df, "DFrame")
        df
    }
