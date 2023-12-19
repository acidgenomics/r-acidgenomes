#' Import Mouse Genome Informatics (MGI) metadata
#'
#' @export
#' @note Updated 2023-12-19.
#'
#' @return `Mgi`.
#'
#' @seealso
#' - https://www.informatics.jax.org/
#'
#' @examples
#' object <- Mgi()
#' print(object)
Mgi <- function() { # nolint
    alert("Importing MGI metadata.")
    url <- pasteUrl(
        "www.informatics.jax.org",
        "downloads",
        "reports",
        "MGI_Gene_Model_Coord.rpt",
        protocol = "https"
    )
    file <- .cacheIt(url)
    lines <- import(con = file, format = "lines")
    cn <- strsplit(lines[[1L]], split = "\t")[[1L]]
    assert(hasLength(cn, n = 15L))
    cn <- sub(pattern = "^[0-9]+\\.\\s", replacement = "", x = cn)
    cn <- camelCase(cn)
    assert(identical(
        x = cn,
        y = c(
            "mgiAccessionId",
            "markerType",
            "markerSymbol",
            "markerName",
            "genomeBuild",
            "entrezGeneId",
            "ncbiGeneChromosome",
            "ncbiGeneStart",
            "ncbiGeneEnd",
            "ncbiGeneStrand",
            "ensemblGeneId",
            "ensemblGeneChromosome",
            "ensemblGeneStart",
            "ensemblGeneEnd",
            "ensemblGeneStrand"
        )
    ))
    cn[cn == "entrezGeneId"] <- "ncbiGeneId"
    cn <- c(cn, "delete")
    lines <- lines[2L:length(lines)]
    con <- textConnection(lines)
    df <- import(
        con = con,
        format = "tsv",
        colnames = cn,
        naStrings = c("NA", "null")
    )
    close(con)
    assert(allAreMatchingFixed(x = df[[1L]], pattern = "MGI:"))
    df <- as(df, "DFrame")
    df[["delete"]] <- NULL
    idCol <- "mgiAccessionId"
    assert(hasNoDuplicates(df[[idCol]]))
    df[[idCol]] <- sub(
        pattern = "^MGI\\:",
        replacement = "",
        x = df[[idCol]]
    )
    df[[idCol]] <- as.integer(df[[idCol]])
    rownames(df) <- df[[idCol]]
    df <- df[order(df[[idCol]]), sort(colnames(df)), drop = FALSE]
    df[["ensemblGeneStrand"]] <- as.factor(df[["ensemblGeneStrand"]])
    df[["ncbiGeneStrand"]] <- as.factor(df[["ncbiGeneStrand"]])
    ## Disabled Rle encoding in 0.7.3 update.
    ## > df <- encode(df)
    metadata(df) <- list(
        "date" = Sys.Date(),
        "organism" = "Mus musculus",
        "packageVersion" = .pkgVersion,
        "url" = url
    )
    new(Class = "Mgi", df)
}
