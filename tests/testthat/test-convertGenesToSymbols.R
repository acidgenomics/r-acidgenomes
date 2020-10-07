context("convertGenesToSymbols")

skip_if_not(hasInternet())

gene2symbol <-
    makeGene2SymbolFromEnsembl(organism = "Homo sapiens", release = 87L)

test_that("character", {
    expect_identical(
        object = convertGenesToSymbols(
            object = c("ENSG00000000003", "ENSG00000000005"),
            gene2symbol = gene2symbol
        ),
        expected = c(
            ENSG00000000003 = "TSPAN6",
            ENSG00000000005 = "TNMD"
        )
    )
})

## Specify organism (to handle FASTA spike-ins (e.g. EGFP).
test_that("FASTA spike-in support", {
    expect_identical(
        object = suppressWarnings({
            convertGenesToSymbols(
                object = c("EGFP", "ENSG00000000003"),
                gene2symbol = gene2symbol
            )
        }),
        expected = c(
            EGFP = "EGFP",
            ENSG00000000003 = "TSPAN6"
        )
    )
})

test_that("Invalid identifiers", {
    expect_warning(
        object = convertGenesToSymbols(
            object = "ENSG00000000000",
            gene2symbol = gene2symbol
        ),
        regexp = "Failed to match genes: ENSG00000000000"
    )
    expect_error(
        object = convertGenesToSymbols(c("ENSG00000000003", NA)),
        regexp = "isCharacter"
    )
    expect_error(
        object = convertGenesToSymbols(c("ENSG00000000003", "")),
        regexp = "isCharacter"
    )
})

test_that("matrix", {
    object <- matrix(
        data = seq_len(4L),
        byrow = TRUE,
        nrow = 2L,
        ncol = 2L,
        dimnames = list(
            c("ENSG00000000003", "ENSG00000000005"),
            c("sample1", "sample2")
        )
    )
    expect_identical(
        object = object %>%
            convertGenesToSymbols(gene2symbol = gene2symbol) %>%
            rownames(),
        expected = c("TSPAN6", "TNMD")
    )
})

test_that("GRanges", {
    object <- gr
    object <- convertGenesToSymbols(object)
    expect_identical(
        object = names(object),
        expected = as.character(Gene2Symbol(gr)[["geneName"]])
    )

})

test_that("SummarizedExperiment", {
    object <- convertGenesToSymbols(rse)
    expect_identical(
        object = rownames(object),
        expected = as.character(mcols(rowRanges(object))[["geneName"]])
    )
})



context("convertSymbolsToGenes")

test_that("SummarizedExperiment", {
    object <- convertGenesToSymbols(rse)
    object <- convertSymbolsToGenes(object)
    expect_identical(
        object = rownames(object),
        expected = as.character(mcols(rowRanges(object))[["geneID"]])
    )
})
