test_that("NcbiGeneHistory", {
    object <- NcbiGeneHistory(organism = "Homo sapiens")
    expect_s4_class(object, "NcbiGeneHistory")
})
