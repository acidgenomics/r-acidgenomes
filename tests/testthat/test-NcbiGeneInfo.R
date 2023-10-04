test_that("NcbiGeneInfo", {
    object <- NcbiGeneInfo(organism = "Homo sapiens")
    expect_s4_class(object, "NcbiGeneInfo")
})
