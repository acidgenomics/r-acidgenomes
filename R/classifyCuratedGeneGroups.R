#' Classify genes into curated HGNC gene groups
#'
#' @export
#' @note Updated 2026-09-24.
#'
#' @details
#' Tags each gene `"riboCyto"` (cytoplasmic ribosomal protein), `"riboMito"`
#' (mitochondrial ribosomal protein), `"hemoglobin"` (hemoglobin subunit), or
#' none of these. Tags are sourced from HGNC's own curated `geneGroupId`
#' assignments, never a symbol regex; a gene not in any curated group gets an
#' empty `character` vector, not an omitted list element.
#'
#' For Mus musculus, human HGNC groups are propagated via a fully
#' identifier-based chain (HGNC `hgncId` -> JAX ortholog `mouseMgiId` -> MGI
#' `ensemblGeneId`), with no gene-symbol matching at any step.
#'
#' This is deliberately independent of `broadClass` (see the internal
#' `.addBroadClass()` in `internal-GenomicRanges.R`): `broadClass` is
#' single-valued and every ribosomal/hemoglobin gene already has a value
#' from it (`"coding"`, `"pseudo"`, etc). Do not fold these tags into
#' `broadClass`.
#'
#' A sibling Python package carries the identical group-ID map and function
#' as a hand-ported twin; a change here must land there too, in the same
#' release.
#'
#' @param ensemblGeneIds `character`.
#' Ensembl gene identifiers to classify.
#'
#' @param organism `character(1)`.
#' Latin organism name. Only `"Homo sapiens"` and `"Mus musculus"` are
#' supported.
#'
#' @param hgnc `Hgnc` or `NULL`.
#' HGNC reference dataset. Downloaded via `Hgnc()` if `NULL`.
#'
#' @param jax `JaxHumanToMouse` or `NULL`.
#' JAX human-to-mouse ortholog dataset. Downloaded via
#' `JaxHumanToMouse(unique = FALSE)` if `NULL`. Ignored for
#' `"Homo sapiens"`.
#'
#' @param mgi `Mgi` or `NULL`.
#' MGI reference dataset. Downloaded via `Mgi()` if `NULL`. Ignored for
#' `"Homo sapiens"`.
#'
#' @return `list`.
#' Named list keyed by `ensemblGeneIds`, each element a `character` vector
#' of curated tags (empty if none).
#'
#' @examples
#' tags <- classifyCuratedGeneGroups(
#'     ensemblGeneIds = "ENSG00000244734",
#'     organism = "Homo sapiens"
#' )
#' print(tags[["ENSG00000244734"]]) # HBB
classifyCuratedGeneGroups <-
    function(
        ensemblGeneIds,
        organism,
        hgnc = NULL,
        jax = NULL,
        mgi = NULL
    ) {
        assert(
            isCharacter(ensemblGeneIds),
            isString(organism)
        )
        if (identical(organism, "Homo sapiens")) {
            tagsByEnsemblId <- .hgncCuratedTagsByEnsemblId(hgnc)
        } else if (identical(organism, "Mus musculus")) {
            tagsByEnsemblId <- .mouseCuratedTagsByEnsemblId(
                hgnc = hgnc,
                jax = jax,
                mgi = mgi
            )
        } else {
            abort(sprintf(
                "Unsupported organism '%s'; only 'Homo sapiens' and 'Mus musculus' are supported.", # nolint
                organism
            ))
        }
        out <- lapply(
            X = ensemblGeneIds,
            FUN = function(id) {
                if (isSubset(id, names(tagsByEnsemblId))) {
                    tagsByEnsemblId[[id]]
                } else {
                    character(0L)
                }
            }
        )
        names(out) <- ensemblGeneIds
        out
    }

## HGNC geneGroupId -> curated tag. Verified live against genenames.org
## 2026-09-24 (hgnc_complete_set.txt): 728 "S ribosomal proteins", 729
## "L ribosomal proteins" (cytoplasmic large/small subunit), 646
## "Mitochondrial ribosomal proteins", 940 "Hemoglobin subunits". The
## "Ribosomal protein S6 kinase family" (1156/1691/3524) is a distinct HGNC
## family (RPS6KA*/RPS6KB*) and is deliberately absent from this map -- do
## not add it.
.curatedGeneGroupIds <- c(
    "728" = "riboCyto",
    "729" = "riboCyto",
    "646" = "riboMito",
    "940" = "hemoglobin"
)

## Resolve one gene's already-split HGNC geneGroupId elements (a plain
## character vector, since Hgnc() returns geneGroupId as a CharacterList)
## to curated tags. NA and non-numeric entries are silently skipped, not
## treated as an error -- most genes carry no gene group at all.
.tagsForGroupIds <- function(groupIds) {
    groupIds <- groupIds[!is.na(groupIds)]
    ids <- suppressWarnings(as.integer(groupIds))
    ids <- ids[!is.na(ids)]
    tags <- unique(unname(.curatedGeneGroupIds[as.character(ids)]))
    tags[!is.na(tags)]
}

## Build human HGNC ID (as character) -> curated tag list from HGNC gene
## groups. Only genes with at least one curated tag are stored.
.tagsByHgncId <- function(hgnc) {
    assert(
        is(hgnc, "Hgnc"),
        isSubset(c("hgncId", "geneGroupId"), colnames(hgnc))
    )
    hgncIds <- as.character(hgnc[["hgncId"]])
    tagsList <- lapply(X = hgnc[["geneGroupId"]], FUN = .tagsForGroupIds)
    keep <- lengths(tagsList) > 0L & !is.na(hgncIds)
    out <- tagsList[keep]
    names(out) <- hgncIds[keep]
    out
}

## Build human Ensembl gene ID -> curated tag list from HGNC gene groups.
.hgncCuratedTagsByEnsemblId <- function(hgnc) {
    if (is.null(hgnc)) {
        hgnc <- Hgnc()
    }
    assert(
        is(hgnc, "Hgnc"),
        isSubset(c("ensemblGeneId", "geneGroupId"), colnames(hgnc))
    )
    ensemblIds <- as.character(hgnc[["ensemblGeneId"]])
    tagsList <- lapply(X = hgnc[["geneGroupId"]], FUN = .tagsForGroupIds)
    keep <- lengths(tagsList) > 0L & !is.na(ensemblIds) & nzchar(ensemblIds)
    out <- tagsList[keep]
    names(out) <- ensemblIds[keep]
    out
}

## Build human HGNC ID -> mouse MGI ID set from the JAX ortholog table.
##
## A human gene can have more than one mouse paralog (e.g. hemoglobin's
## Hbb-bs/Hbb-bt/Hbb-bh2/Hbb-b1/Hbb-b2 expansion) -- every mouse MGI ID
## reachable from a given human HGNC ID is collected, not just the
## first/last one seen. `JaxHumanToMouse()` never has this repo's Python
## sibling's mouse_mgi_id column-collision bug: it explicitly drops the
## meaningless opposite-species copy of each shared field before merging
## (see its own `hs[["mouseMgiId"]] <- NULL` / `mm[["hgncId"]] <- NULL`).
.mouseMgiIdsByHumanHgncId <- function(jax) {
    assert(
        is(jax, "JaxHumanToMouse"),
        isSubset(c("humanHgncId", "mouseMgiId"), colnames(jax))
    )
    humanHgncIds <- as.character(jax[["humanHgncId"]])
    mouseMgiIds <- jax[["mouseMgiId"]]
    keep <- !is.na(humanHgncIds) & !is.na(mouseMgiIds)
    split(x = mouseMgiIds[keep], f = humanHgncIds[keep])
}

## Build mouse MGI ID (as character) -> Ensembl gene ID from the MGI
## gene-model report.
.mouseEnsemblIdByMgiId <- function(mgi) {
    assert(
        is(mgi, "Mgi"),
        isSubset(c("mgiAccessionId", "ensemblGeneId"), colnames(mgi))
    )
    mgiIds <- as.character(mgi[["mgiAccessionId"]])
    ensemblIds <- as.character(mgi[["ensemblGeneId"]])
    keep <- !is.na(mgiIds) & !is.na(ensemblIds) & nzchar(ensemblIds)
    out <- as.list(ensemblIds[keep])
    names(out) <- mgiIds[keep]
    out
}

## Propagate human HGNC curated tags to mouse Ensembl gene IDs.
##
## Fully identifier-based: HGNC `hgncId` -> JAX `humanHgncId` /
## `mouseMgiId` -> MGI `mgiAccessionId` / `ensemblGeneId`. No gene-symbol
## matching at any step, unlike a naive mouse-symbol regex (which silently
## returns zero matches -- hemoglobin symbols are hyphenated in mouse,
## e.g. `"Hba-a1"`, and share no substring with any human-derived pattern).
.mouseCuratedTagsByEnsemblId <- function(hgnc, jax, mgi) {
    if (is.null(hgnc)) {
        hgnc <- Hgnc()
    }
    if (is.null(jax)) {
        jax <- JaxHumanToMouse(unique = FALSE)
    }
    if (is.null(mgi)) {
        mgi <- Mgi()
    }
    tagsByHgncId <- .tagsByHgncId(hgnc)
    mgiIdsByHgncId <- .mouseMgiIdsByHumanHgncId(jax)
    ensemblIdByMgiId <- .mouseEnsemblIdByMgiId(mgi)

    out <- list()
    for (hgncId in names(tagsByHgncId)) {
        mgiIds <- as.character(mgiIdsByHgncId[[hgncId]])
        for (mgiId in mgiIds) {
            ensemblId <- ensemblIdByMgiId[[mgiId]]
            if (!is.null(ensemblId)) {
                out[[ensemblId]] <- tagsByHgncId[[hgncId]]
            }
        }
    }
    out
}
