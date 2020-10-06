#' Map protein identifiers to genes
#'
#' @name makeProtein2Gene
#' @note Updated 2020-10-06.
#'
#' @inheritParams AcidRoxygen::params
#' @param ids `character`.
#'   Ensembl protein identifiers.
#'   Human proteins are prefixed with "ENSP", for example.
#'
#' @return `Protein2Gene`.
#'
#' @examples
#' ids <- c("ENSP00000238714", "ENSP00000338157")
#' makeProtein2GeneFromEnsembl(ids)
NULL



#' @rdname makeProtein2Gene
#' @export
makeProtein2GeneFromEnsembl <- function(
    ids,
    organism = NULL,
    genomeBuild = NULL,
    release = NULL
) {
    assert(
        isCharacter(ids),
        hasNoDuplicates(ids)
    )
    cli_alert("Making {.var Protein2Gene} from Ensembl.")
    if (is.null(organism)) {
        organism <- detectOrganism(ids)
    }
    edb <- getEnsDb(
        organism = organism,
        genomeBuild = genomeBuild,
        release = release
    )
    df <- select(
        x = edb,
        keys = ids,
        keytype = "PROTEINID",
        columns = c("GENEID", "GENENAME")
    )
    df <- as(df, "DataFrame")
    colnames(df) <- camelCase(colnames(df))
    colnames(df) <- gsub("id$", "ID", colnames(df))
    colnames(df) <- gsub("name$", "Name", colnames(df))
    if (!areSetEqual(ids, unique(df[["proteinID"]]))) {
        stop(sprintf(
            "Match failure: %s.",
            toString(
                sort(setdiff(ids, unique(df[["proteinID"]]))),
                width = 500L
            )
        ))
    }
    metadata(df) <- .getEnsDbMetadata(edb)
    new(Class = "Protein2Gene", df)
}
