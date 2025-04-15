startsWith <- IRanges::startsWith

hgnc <- Hgnc()
i <- !is.na(hgnc[["ensemblGeneId"]])
hgnc <- hgnc[i, ]
hgncGeneNames <- unlist(
    x = list(
        hgnc[["geneName"]],
        unlist(hgnc[["aliasSymbol"]]),
        unlist(hgnc[["prevSymbol"]])
    ),
    recursive = TRUE,
    use.names = FALSE
)
hgncGeneNames <- sort(unique(hgncGeneNames))

ncbi <- NcbiGeneInfo(organism = "Homo sapiens")
i <- any(startsWith(x = ncbi[["dbXrefs"]], prefix = "Ensembl:"))
## NOTE Consider removing NCBI genes that multi-map to multiple Ensembl IDs.
ncbi <- ncbi[i, ]
ncbiGeneNames <- unlist(
    x = list(
        ncbi[["geneName"]]
        ## > unlist(ncbi[["geneSynonyms"]])
    ),
    recursive = TRUE,
    use.names = FALSE
)
ncbiGeneNames <- sort(unique(ncbiGeneNames))

test_that("HGNC gene names", {
    expect_no_error(
        object = mapGeneNamesToEnsembl(
            genes = hgncGeneNames,
            organism = "Homo sapiens",
            hgnc = hgnc
        )
    )
    expect_error(
        object = mapGeneNamesToEnsembl(
            genes = hgncGeneNames,
            organism = "Homo sapiens",
            ncbi = ncbi
        ),
        regexp = "mapping failure"
    )
})

## FIXME Rework this to use direct mappings instead.
## FIXME: "AARS1P1", "AARSD1P1", "ABCA15P", "ABCB10P3"
## 106480683 to ENSG00000249038
test_that("NCBI gene names", {
    out <- mapGeneNamesToEnsembl(
        genes = ncbiGeneNames,
        organism = "Homo sapiens",
        ncbi = ncbi
    )
})




out <- mapGeneNamesToEnsembl(
    genes = ncbiGeneNames,
    organism = "Homo sapiens",
    hgnc = hgnc
)
## ! 176965 mapping failures: "'C-K-RAS",
## "(FM-3)", "(IV)-44", "(ppGpp)ase",
## "0610011B16Rik", "0610037N12Rik",
## "0710008D09Rik", "1-12P", "1-14P", "1-17P",
## "1-67P", "1-68P", "1-AGPAT 3", "1-AGPAT 6",
## "1-AGPAT1", "1-AGPAT2", "1-AGPAT4",
## "10-FTHFDH", "104p", "105A"....




## HGNC Currently has "HDGFRP3" incorrectly labeled as "Hdgfrp3" (note case).
# ! 1 mapping failure: "HDGFRP3"

test_that("Case sensitivity", {
    expect_identical(
        object = mapGeneNamesToEnsembl(
            genes = c("HDGFL3", "Hdgfrp3"),
            organism = "Homo sapiens",
            ignoreCase = TRUE,
            hgnc = hgnc
        ),
        expected = rep("ENSG00000166503", times = 2L)
    )
    expect_error(
        object = mapGeneNamesToEnsembl(
            genes = c("HDGFL3", "HDGFRP3"),
            organism = "Homo sapiens",
            ignoreCase = FALSE,
            hgnc = hgnc
        ),
        regexp = "HDGFRP3"
    )
})




## Works.
mapGeneNamesToEnsembl(
    genes = c("HDGFL3", "HDGFRP3"),
    organism = "Homo sapiens",
    ncbi = ncbi
)
