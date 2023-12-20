test_that("Homo sapiens : exons", {
    object <- makeGRangesFromEnsembl(
        organism = "Homo sapiens",
        level = "exons",
        release = 110L,
        ignoreVersion = TRUE
    )
    expect_s4_class(object, "EnsemblExons")
    expect_length(object, 735689L)
    ## ensembldb doesn't currently support versioned exon identifiers.
    expect_error(
        object = makeGRangesFromEnsembl(
            organism = "Homo sapiens",
            level = "exons",
            ignoreVersion = FALSE
        ),
        regexp = "ignoreVersion"
    )
})

test_that("Mus musculus : exons", {
    object <- makeGRangesFromEnsembl(
        organism = "Mus musculus",
        level = "exons",
        release = 110L,
        ignoreVersion = TRUE
    )
    expect_s4_class(object, "EnsemblExons")
    expect_length(object, 458023L)
})
