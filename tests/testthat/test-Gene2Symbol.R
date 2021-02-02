context("Gene2Symbol")

formats <- eval(formals(`Gene2Symbol,GRanges`)[["format"]])
test_that("Gene2Symbol", {
    for (format in formats) {
        object <- GRanges
        object <- Gene2Symbol(object, format = format)
        expect_s4_class(object, "Gene2Symbol")
        expect_identical(colnames(object), c("geneId", "geneName"))
    }
})
rm(formats)

test_that("summary", {
    object <- GRanges
    x <- Gene2Symbol(object)
    output <- capture.output(summary(x))
    expect_identical(
        head(output, n = 1L),
        paste("genes:", nrow(object))
    )
})
