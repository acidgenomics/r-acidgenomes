#' Map input to human gene orthologs
#'
#' @export
#' @note Updated 2021-08-03.
#'
#' @details
#' Genes with identifier versions (e.g. "ENSMUSG00000000001.5") are not
#' currently supported.
#'
#' @inheritParams AcidRoxygen::params
#'
#' @return `DataFrame`.
#' Data frame containing mapping columns:
#'
#' - `geneId`
#' - `geneName`
#' - `humanGeneId`
#' - `humanGeneName`
#'
#' @seealso
#' - `biomaRt::listEnsemblArchives()`.
#' - `biomaRt::listMarts()`.
#' - `biomaRt::useMart()`.
#'
#' @examples
#' ## Expecting ENSMUSG00000000001, ENSMUSG00000000003 to not match here.
#' genes <- c(
#'     "ENSMUSG00000000001", "ENSMUSG00000000003",
#'     "ENSMUSG00000000028", "ENSMUSG00000000031",
#'     "ENSMUSG00000000037", "ENSMUSG00000000049"
#' )
#' ## Protect against Ensembl timeouts causing build checks to fail.
#' if (goalie::hasInternet("https://ensembl.org")) {
#'     try({
#'         x <- mapHumanOrthologs(genes = genes, ensemblRelease = 87L)
#'         print(x)
#'     })
#' }
mapHumanOrthologs <-
    function(genes,
             organism = NULL,
             ensemblRelease = NULL) {
        pkgs <- .packages()
        requireNamespaces("biomaRt")
        assert(
            isCharacter(genes),
            isOrganism(organism, nullOK = TRUE),
            isInt(ensemblRelease, nullOK = TRUE)
        )
        if (is.null(organism)) {
            organism <- detectOrganism(genes)
        }
        ## Don't allow the user to pass in human genes.
        assert(!identical(organism, "Homo sapiens"))
        ## Match the Ensembl release to the archive host name.
        ## e.g. Ensembl 99: http://jan2020.archive.ensembl.org
        host <- mapEnsemblReleaseToURL(ensemblRelease)
        ## e.g. "mmusculus_gene_ensembl".
        dataset <- paste0(
            tolower(sub(
                pattern = "^([A-Z])([a-z]+)\\s([a-z]+)$",
                replacement = "\\1\\3",
                x = organism
            )),
            "_gene_ensembl"
        )
        alert(sprintf(
            fmt = paste(
                "Matching orthologs against {.var %s} ({.url %s}) with",
                "{.pkg %s} %s."
            ),
            dataset, host, "bioMart",
            as.character(packageVersion("biomaRt"))
        ))
        ## Can use "ENSEMBL_MART_ENSEMBL" instead of "ensembl" here.
        mart <- tryCatch(
            expr = biomaRt::useMart(
                biomart = "ensembl",
                dataset = dataset,
                host = host,
                verbose = FALSE
            ),
            error = function(e) {
                abort(sprintf(
                    "{.pkg %s}::{.fun %s} failure.",
                    "biomaRt", "useMart"
                ))
            }
        )
        map <- tryCatch(
            expr = biomaRt::select(
                x = mart,
                keys = genes,
                keytype = "ensembl_gene_id",
                columns = c(
                    "ensembl_gene_id",
                    "hsapiens_homolog_ensembl_gene"
                )
            ),
            error = function(e) {
                abort(sprintf(
                    "{.pkg %s}::{.fun %s} failure.",
                    "biomaRt", "select"
                ))
            }
        )
        if (!hasRows(map)) {
            abort("Failed to map any genes.")
        }
        map <- as(map, "DataFrame")
        colnames(map) <- c("geneId", "humanGeneId")
        map <- sanitizeNA(map)
        keep <- complete.cases(map)
        if (!all(keep)) {
            n <- sum(!keep)
            failures <- sort(unique(map[which(!keep), "geneId"]))
            alertWarning(sprintf(
                "%d match %s: %s",
                n,
                ngettext(
                    n = n,
                    msg1 = "failure",
                    msg2 = "failures"
                ),
                toInlineString(failures, n = 5L)
            ))
            map <- map[keep, , drop = FALSE]
        }
        map <- map[order(map), , drop = FALSE]
        keep <- !duplicated(map[["geneId"]])
        if (!all(keep)) {
            dupes <- unique(map[["geneId"]][!keep])
            alertWarning(sprintf(
                paste(
                    "%d gene %s multimap to multiple human genes.",
                    "Selecting first match by oldest Ensembl identifier."
                ),
                length(dupes),
                ngettext(
                    n = length(dupes),
                    msg1 = "identifier",
                    msg2 = "identifiers"
                )
            ))
            map <- map[keep, , drop = FALSE]
        }
        keep <- !duplicated(map[["humanGeneId"]])
        if (!all(keep)) {
            dupes <- unique(map[["humanGeneId"]][!keep])
            alertWarning(sprintf(
                paste(
                    "%d duplicate human gene %s detected.",
                    "Selecting first match."
                ),
                length(dupes),
                ngettext(
                    n = length(dupes),
                    msg1 = "identifier",
                    msg2 = "identifiers"
                )
            ))
            map <- map[keep, , drop = FALSE]
        }
        assert(
            hasRows(map),
            hasNoDuplicates(map[["geneId"]]),
            hasNoDuplicates(map[["humanGeneId"]]),
            isTRUE(all(complete.cases(map)))
        )
        alertInfo(sprintf(
            "%d gene %s mapped 1:1 from {.emph %s} to {.emph %s}.",
            nrow(map),
            ngettext(
                n = nrow(map),
                msg1 = "identifier",
                msg2 = "identifiers"
            ),
            organism, "Homo sapiens"
        ))
        alert(sprintf("Getting {.emph %s} gene symbols.", organism))
        g2s <- makeGene2SymbolFromEnsembl(
            organism = organism,
            release = ensemblRelease,
            ignoreVersion = TRUE,
            format = "unmodified"
        )
        g2s <- as(g2s, "DataFrame")
        assert(identical(colnames(g2s), c("geneId", "geneName")))
        alert(sprintf("Getting {.emph %s} gene symbols.", "Homo sapiens"))
        g2sHuman <- makeGene2SymbolFromEnsembl(
            organism = "Homo sapiens",
            release = ensemblRelease,
            ignoreVersion = TRUE,
            format = "unmodified"
        )
        g2sHuman <- as(g2sHuman, "DataFrame")
        colnames(g2sHuman) <- c("humanGeneId", "humanGeneName")
        assert(all(complete.cases(g2sHuman)))
        ## Return.
        out <- map
        out <- leftJoin(out, g2s, by = "geneId")
        out <- leftJoin(out, g2sHuman, by = "humanGeneId")
        rownames(out) <- out[["geneId"]]
        out <- out[, sort(colnames(out))]
        assert(
            !all(is.na(out[["geneName"]])),
            !all(is.na(out[["humanGeneId"]]))
        )
        forceDetach(keep = pkgs)
        out
    }
