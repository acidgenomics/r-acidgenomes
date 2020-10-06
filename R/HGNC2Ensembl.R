#' @inherit HGNC2Ensembl-class title description return
#' @note Updated 2020-10-03.
#' @export
#' @examples
#' object <- HGNC2Ensembl()
#' print(object)
HGNC2Ensembl <-  # nolint
    function() {
        hgnc <- HGNC()
        cli_alert("Mapping HGNC identifiers to Ensembl.")
        df <- as(hgnc, "DataFrame")
        cols <- c("hgncID", "ensemblGeneID")
        assert(isSubset(cols, colnames(df)))
        df <- df[, cols]
        colnames(df)[colnames(df) == "ensemblGeneID"] <- "ensembl"
        colnames(df)[colnames(df) == "hgncID"] <- "hgnc"
        df <- df[complete.cases(df), ]
        metadata(df) <- metadata(hgnc)
        new(Class = "HGNC2Ensembl", df)
    }
