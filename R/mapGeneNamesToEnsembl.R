#' Map gene names to Ensembl identifiers
#'
#' @export
#' @note Updated 2025-04-15.
#'
#' @details Internally matches using `mapGeneNamesToHgnc` (*Homo sapiens* only)
#' or `mapGeneNamesToNcbi` (all other organisms), so we can support gene synonym
#' matching.
#'
#' @inheritParams AcidRoxygen::params
#'
#' @param genes
#' Gene names (e.g. `"TUT4"`).
#'
#' @param hgnc `Hgnc` object.
#' Supported for *Homo sapiens* genome only.
#' Snapshot of HGNC annotations.
#' Passes to `mapGeneNamesToHgnc` internally.
#'
#' @param ncbi `NcbiGeneInfo` object.
#' Snapshot of NCBI annotations.
#' Passes to `mapGeneNamesToNcbi` internally.
#'
#' @examples
#' ## Homo sapiens.
#' x <- mapGeneNamesToEnsembl(
#'     genes = c("TUT4", "ZCCHC11", "TENT3A"),
#'     organism = "Homo sapiens"
#' )
#' print(x)
#'
#' ## Mus musculus
#' x <- mapGeneNamesToEnsembl(
#'     genes = c("Nfe2l2", "Nrf2"),
#'     organism = "Mus musculus"
#' )
#' print(x)
mapGeneNamesToEnsembl <-
    function(genes,
             organism,
             ignoreCase = FALSE,
             hgnc = NULL,
             ncbi = NULL) {
        assert(
            isCharacter(genes),
            isOrganism(organism),
            isFlag(ignoreCase),
            is(ncbi, "NcbiGeneInfo") || is.null(ncbi),
            is(hgnc, "Hgnc") || is.null(hgnc)
        )
        ## Default to HGNC over NCBI for Homo sapiens.
        if (
            identical(organism, "Homo sapiens") &&
            is.null(hgnc) &&
            is.null(ncbi)
        ) {
            hgnc <- Hgnc()
        }
        if (is(hgnc, "Hgnc")) {
            assert(
                identical(organism, "Homo sapiens"),
                is.null(ncbi)
            )
            keep <- !is.na(hgnc[["ensemblGeneId"]])
            hgnc <- hgnc[keep, ]
            ids <- mapGeneNamesToHgnc(
                genes = genes,
                ignoreCase = ignoreCase,
                hgnc = hgnc
            )
            map <- as(hgnc, "DFrame")
            map <- map[, c("hgncId", "ensemblGeneId")]
        } else {
            assert(is.null(hgnc))
            keep <- any(startsWith(x = ncbi[["dbXrefs"]], prefix = "Ensembl:"))
            ncbi <- ncbi[keep, ]
            ensGene <- ncbi[["dbXrefs"]]
            ensGene <- ensGene[startsWith(x = ensGene, prefix = "Ensembl:")]
            keep <- lengths(ensGene) == 1L
            ensGene <- ensGene[keep]
            ncbi <- ncbi[keep, ]
            ids <- mapGeneNamesToNcbi(
                genes = genes,
                organism = organism,
                ignoreCase = ignoreCase,
                ncbi = ncbi
            )
            map <- as(ncbi, "DFrame")
            ensGene <- sub(pattern = "^Ensembl:", replacement = "", x = ensGene)
            ensGene <- unlist(ensGene, recursive = FALSE, use.names = FALSE)
            map[["ensemblGeneId"]] <- ensGene
            map <- map[, c("geneId", "ensemblGeneId")]
        }
        idx <- match(x = ids, table = map[[1L]])
        if (anyNA(idx)) {
            fail <- genes[is.na(idx)]
            abort(sprintf(
                "%d mapping %s: %s.",
                length(fail),
                ngettext(
                    n = length(fail),
                    msg1 = "failure",
                    msg2 = "failures"
                ),
                toInlineString(x = fail, n = 20L)
            ))
        }
        out <- map[idx, 2L, drop = TRUE]
        out
    }
