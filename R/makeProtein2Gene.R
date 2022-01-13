#' Map protein identifiers to genes
#'
#' @name makeProtein2Gene
#' @note Updated 2022-01-12.
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
#' p2g <- makeProtein2GeneFromEnsembl(ids)
#' print(p2g)
NULL



#' @rdname makeProtein2Gene
#' @export
makeProtein2GeneFromEnsembl <- function(
    ids,
    organism = NULL,
    genomeBuild = NULL,
    release = NULL
) {
    pkgs <- .packages()
    requireNamespaces("ensembldb")
    assert(
        isCharacter(ids),
        hasNoDuplicates(ids)
    )
    alert(sprintf("Making {.cls %s} from Ensembl.", "Protein2Gene"))
    if (is.null(organism)) {
        organism <- detectOrganism(ids)
    }
    edb <- .getEnsDb(
        organism = organism,
        genomeBuild = genomeBuild,
        release = release
    )
    ## The `select()` generic is defined in AnnotationDbi.
    df <- ensembldb::select(
        x = edb,
        keys = ids,
        keytype = "PROTEINID",
        columns = c("GENEID", "GENENAME")
    )
    df <- as(df, "DataFrame")
    colnames(df) <- tolower(colnames(df))
    colnames(df) <- gsub("id$", "Id", colnames(df))
    colnames(df) <- gsub("name$", "Name", colnames(df))
    if (!areSetEqual(ids, unique(df[["proteinId"]]))) {
        abort(sprintf(
            "Match failure: %s.",
            toInlineString(
                sort(setdiff(ids, unique(df[["proteinId"]]))),
                n = 10L
            )
        ))
    }
    metadata(df) <- .getEnsDbMetadata(edb)
    forceDetach(keep = pkgs)
    new(Class = "Protein2Gene", df)
}
