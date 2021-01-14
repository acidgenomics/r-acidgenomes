#' @inherit HGNC2Ensembl-class title description return
#' @note Updated 2021-01-06.
#' @export
#' @examples
#' object <- HGNC2Ensembl()
#' print(object)
HGNC2Ensembl <-  # nolint
    function() {
        hgnc <- HGNC()
        alert("Mapping HGNC identifiers to Ensembl.")
        df <- as(hgnc, "DataFrame")
        cols <- c("hgncId", "ensemblGeneId")
        assert(isSubset(cols, colnames(df)))
        df <- df[, cols]
        colnames(df)[colnames(df) == "ensemblGeneId"] <- "ensembl"
        colnames(df)[colnames(df) == "hgncId"] <- "hgnc"
        df <- df[complete.cases(df), ]
        metadata(df) <- metadata(hgnc)
        new(Class = "HGNC2Ensembl", df)
    }
