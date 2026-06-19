test_that("Homo sapiens : GENCODE", {
    gr <- makeGRangesFromGff(
        file = gffs[["gencode_grch38_gtf"]],
        ignoreVersion = FALSE,
        extraMcols = FALSE
    )
    g2s <- GeneToSymbol(object = gr, format = "1:1")
    expect_s4_class(g2s, "GeneToSymbol")
    expect_identical(nrow(g2s), 61228L)
})

test_that("Homo sapiens : RefSeq", {
    grl <- makeGRangesFromGff(
        file = gffs[["refseq_grch38_gtf"]],
        ignoreVersion = FALSE,
        extraMcols = FALSE
    )
    g2s <- GeneToSymbol(object = grl, format = "1:1")
    expect_s4_class(g2s, "GeneToSymbol")
    expect_identical(nrow(g2s), 59638L)
})
