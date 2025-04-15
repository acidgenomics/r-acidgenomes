test_that("Homo sapiens", {
    df <- importEnsemblToNcbiGeneMap(
        organism = "Homo sapiens",
        genomeBuild = "GRCh38",
        release = 113L,
        uniqueOnly = TRUE
    )
    expect_s4_class(df, "DFrame")
    expect_identical(nrow(df), 25226L)
    expect_identical(df[1L, "ensemblGeneId"], "ENSG00000000003")
    expect_identical(df[1L, "ncbiGeneId"], 7105L)
})

test_that("Mus musculus", {
    df <- importEnsemblToNcbiGeneMap(
        organism = "Mus musculus",
        genomeBuild = "GRCm39",
        release = 113L,
        uniqueOnly = TRUE
    )
    expect_s4_class(df, "DFrame")
    expect_identical(nrow(df), 27934L)
    expect_identical(df[1L, "ensemblGeneId"], "ENSMUSG00000000001")
    expect_identical(df[1L, "ncbiGeneId"], 14679L)
})
