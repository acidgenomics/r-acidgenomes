#' Import human-to-mouse gene mappings
#'
#' @export
#' @note Updated 2023-09-26.
#'
#' @return `HumanToMouse`.
#'
#' @seealso
#' - https://www.informatics.jax.org/
#' - https://www.biostars.org/p/9567892/
#'
#' @examples
#' h2m <- HumanToMouse()
#' print(h2m)
HumanToMouse <- function() {
    url <- pasteURL(
        "www.informatics.jax.org",
        "downloads",
        "reports",
        "HOM_MouseHumanSequence.rpt",
        protocol = "https"
    )
    file <- .cacheIt(url)
    df <- import(file, format = "tsv")
    df <- as(df, "DFrame")
    colnames(df) <- camelCase(colnames(df))
    cols <- c(
        "dbClassKey", "entrezGeneId", "hgncId", "ncbiTaxonId", "mouseMgiId",
        "omimGeneId", "symbol"
    )
    assert(isSubset(cols, colnames(df)))
    df <- df[, cols]
    df[["hgncId"]] <- sub(
        pattern = "^HGNC:",
        replacement = "",
        x = df[["hgncId"]]
    )
    df[["hgncId"]] <- as.integer(df[["hgncId"]])
    df[["mouseMgiId"]] <- sub(
        pattern = "^MGI:",
        replacement = "",
        x = df[["mouseMgiId"]]
    )
    df[["mouseMgiId"]] <- as.integer(df[["mouseMgiId"]])
    df[["omimGeneId"]] <- sub(
        pattern = "^OMIM:",
        replacement = "",
        x = df[["omimGeneId"]]
    )
    df[["omimGeneId"]] <- as.integer(df[["omimGeneId"]])
    spl <- split(df, f = df[["ncbiTaxonId"]])
    assert(is(spl, "SplitDFrameList"))
    hs <- spl[["9606"]]
    hs[["ncbiTaxonId"]] <- NULL
    hs[["mouseMgiId"]] <- NULL
    colnames(hs)[colnames(hs) == "entrezGeneId"] <- "humanNcbiGeneId"
    colnames(hs)[colnames(hs) == "hgncId"] <- "humanHgncId"
    colnames(hs)[colnames(hs) == "omimGeneId"] <- "humanOmimGeneId"
    colnames(hs)[colnames(hs) == "symbol"] <- "humanGeneName"
    mm <- spl[["10090"]]
    mm[["hgncId"]] <- NULL
    mm[["ncbiTaxonId"]] <- NULL
    mm[["omimGeneId"]] <- NULL
    colnames(mm)[colnames(mm) == "entrezGeneId"] <- "mouseNcbiGeneId"
    colnames(mm)[colnames(mm) == "symbol"] <- "mouseGeneName"
    df <- merge(x = hs, y = mm, by = "dbClassKey", all.x = TRUE, all.y = TRUE)
    assert(is(df, "DFrame"))
    keep <- complete.cases(df[, c("humanGeneName", "mouseGeneName")])
    df <- df[keep, sort(colnames(df))]
    idx <- order(df[["humanGeneName"]], df[["mouseGeneName"]])
    df <- df[idx, ]
    df
}
