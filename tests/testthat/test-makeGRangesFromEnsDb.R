## Conditionally test if optional EnsDb.Hsapiens.v75 package is installed.
pkg <- "EnsDb.Hsapiens.v75"
skip_if_not_installed(pkg)

test_that("Genes", {
    x <- makeGRangesFromEnsDb(object = pkg, level = "genes")
    expect_s4_class(x, "GenomicRanges")
})

## Note that transcript versions aren't saved in this object.
test_that("Transcripts", {
    x <- makeGRangesFromEnsDb(
        object = pkg,
        level = "transcripts",
        ignoreVersion = TRUE
    )
    expect_s4_class(x, "GenomicRanges")
})

ah <- AnnotationHub::AnnotationHub(ask = FALSE)
edb <- ah[["AH73881"]]

test_that("Genes", {
    x <- makeGRangesFromEnsDb(object = edb, level = "genes")
    expect_s4_class(x, "GenomicRanges")
})

test_that("Transcripts", {
    x <- makeGRangesFromEnsDb(
        object = edb,
        level = "transcripts",
        ignoreVersion = FALSE
    )
    expect_s4_class(x, "GenomicRanges")
})
