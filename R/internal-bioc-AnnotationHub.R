#' Connect to AnnotationHub
#'
#' @note Updated 2020-09-24.
#' @noRd
#'
#' @details
#' On a fresh install this will print a txProgressBar to the console. We're
#' using [utils::capture.output()] here to suppress the console output, since
#' it's not very informative and can cluster R Markdown reports.
.annotationHub <- function() {
    userAttached <- .packages()
    invisible(capture.output({
        ah <- suppressMessages(AnnotationHub())
    }))
    assert(is(ah, "AnnotationHub"))
    forceDetach(keep = userAttached)
    ah
}



#' Get Ensembl/Entrez mappings from NCBI OrgDb via AnnotationHub
#'
#' @note Updated 2020-10-01.
#' @noRd
.getEnsembl2EntrezFromOrgDb <- function(
    keys,
    keytype,
    columns,
    organism,
    strict = TRUE
) {
    assert(
        isCharacter(keys),
        hasNoDuplicates(keys),
        isString(keytype),
        isCharacter(columns),
        isOrganism(organism),
        isFlag(strict)
    )
    alert(sprintf(
        "Matching identifiers using NCBI {.var %s} via {.pkg %s} %s.",
        "OrgDb",
        "AnnotationHub",
        packageVersion("AnnotationHub")
    ))
    userAttached <- .packages()
    ah <- .annotationHub()
    ahs <- query(ah, pattern = c(organism, "NCBI", "OrgDb"))
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
    suppressMessages({
        df <- select(
            x = orgdb,
            keys = keys,
            keytype = keytype,
            columns = columns
        )
    })
    df <- as(df, "DataFrame")
    if (isTRUE(strict)) {
        df <- df[complete.cases(df), ]
        if (!areSetEqual(keys, unique(df[[keytype]]))) {
            setdiff <- setdiff(keys, unique(df[[keytype]]))
            stop(sprintf(
                "Match failure: %s.",
                toString(setdiff, width = 200L)
            ))
        }
    }
    colnames(df)[colnames(df) == "ENSEMBL"] <- "ensembl"
    colnames(df)[colnames(df) == "ENTREZID"] <- "entrez"
    df[["entrez"]] <- as.integer(df[["entrez"]])
    forceDetach(keep = userAttached)
    df
}
