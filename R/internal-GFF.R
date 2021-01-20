## Updated 2021-01-20.
.gffMetadataForTxDb <- function(file) {
    df <- getGFFMetadata(file, nMax = 2000L)
    if (is.null(db)) return(NULL)
    list(
        "genomeBuild" = .gffGenomeBuild(df),
        "source" = .gffSource(df)
    )
}



## Updated 2021-01-20.
.gffGenomeBuild <- function(df) {
    assert(is(df, "DataFrame"))
    ## GENCODE files have a description key that contains the genome build.
    if (isTRUE("description" %in% df[["key"]])) {
        string <- df[df[["key"]] == "description", "value", drop = TRUE]
        x <- str_match(
            string = string,
            pattern = "genome \\(([^\\)]+)\\)"
        )[1L, 2L]
        if (isString(x)) return(x)
    }
    ## Otherwise we can parse for standard "genome-build" key, which is
    ## supported by Ensembl and RefSeq.
    .getValue <- function(key) {
        x <- df[match(x = key, table = df[["key"]]), "value", drop = TRUE]
        if (is.na(x)) return(NULL)
        x
    }
    x <- .getValue("genome-build")
    if (isString(x)) return(x)
    NULL
}



## Updated 2021-01-20.
.gffSource <- function(df) {
    assert(is(df, "DataFrame"))
    annoSource <-
        df[df[["key"]] == "annotation-source", "value", drop = TRUE]
    provider <-
        df[df[["key"]] == "provider", "value", drop = TRUE]
    if (provider == "GENCODE") {
        return("GENCODE")
    } else if (isTRUE(grepl(pattern = "^NCBI", x = annoSource))) {
        return("RefSeq")
    } else if (isSubset(
        x = c(
            "genebuild-last-updated",
            "genome-build",
            "genome-build-accession",
            "genome-date",
            "genome-version"
        ),
        y = df[["key"]]
    )) {
        return("Ensembl")
    }
    NULL
}
