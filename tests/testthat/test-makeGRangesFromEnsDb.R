context("makeGRangesFromEnsDb : GRCh37 (hg19)")

## Conditionally test if optional EnsDb.Hsapiens.v75 package is installed.
pkg <- "EnsDb.Hsapiens.v75"
skip_if_not_installed(pkg)

test_that("Genes", {
    x <- makeGRangesFromEnsDb(object = pkg, level = "genes")
    expect_s4_class(x, "GRanges")
})

## Note that transcript versions aren't saved in this object.
test_that("Transcripts", {
    x <- makeGRangesFromEnsDb(
        object = pkg,
        level = "transcripts",
        ignoreVersion = TRUE
    )
    expect_s4_class(x, "GRanges")
})



context("makeGRangesFromEnsDb : AnnotationHub")

ah <- AnnotationHub::AnnotationHub()
edb <- ah[["AH73881"]]

test_that("Genes", {
    x <- makeGRangesFromEnsDb(object = edb, level = "genes")
    expect_s4_class(x, "GRanges")
})

test_that("Transcripts", {
    x <- makeGRangesFromEnsDb(
        object = edb,
        level = "transcripts",
        ignoreVersion = FALSE
    )
    expect_s4_class(x, "GRanges")
})
