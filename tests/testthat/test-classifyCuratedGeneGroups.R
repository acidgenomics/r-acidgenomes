test_that("Homo sapiens: named gene memberships", {
    tags <- classifyCuratedGeneGroups(
        ensemblGeneIds = c(
            "ENSG00000244734", # HBB
            "ENSG00000225323", # HBAP1, a pseudogene -- HGNC's own group
            #  includes it
            "ENSG00000198034", # RPS4X
            "ENSG00000169288", # MRPL1
            "ENSG00000117676", # RPS6KA1: S6 kinase family, must NOT be
            #  tagged ribosomal
            "ENSG00000000000" # not a real Ensembl ID
        ),
        organism = "Homo sapiens"
    )
    expect_identical(tags[["ENSG00000244734"]], "hemoglobin")
    expect_identical(tags[["ENSG00000225323"]], "hemoglobin")
    expect_identical(tags[["ENSG00000198034"]], "riboCyto")
    expect_identical(tags[["ENSG00000169288"]], "riboMito")
    expect_identical(tags[["ENSG00000117676"]], character(0L))
    expect_identical(tags[["ENSG00000000000"]], character(0L))
})

test_that("Mus musculus: fully identifier-based propagation", {
    ## Confirms the mouse path never falls back to gene-symbol matching --
    ## a naive human-derived regex would return zero hemoglobin genes for
    ## mouse, since mouse hemoglobin symbols are hyphenated (e.g.
    ## "Hba-a1") and share no substring with any human pattern.
    tags <- classifyCuratedGeneGroups(
        ensemblGeneIds = "ENSMUSG00000052305", # Hbb-bs
        organism = "Mus musculus"
    )
    expect_identical(tags[["ENSMUSG00000052305"]], "hemoglobin")
})

test_that("unsupported organism raises", {
    expect_error(
        object = classifyCuratedGeneGroups(
            ensemblGeneIds = "ENSG001",
            organism = "Danio rerio"
        ),
        regexp = "Unsupported organism"
    )
})
