context("convertUCSCBuildToEnsembl")

test_that("convertUCSCBuildToEnsembl", {
    mapply(
        object = list(
            "hg19",
            "hg38",
            "mm10"
        ),
        expected = list(
            "GRCh37",
            "GRCh38",
            "GRCm38"
        ),
        FUN = function(object, expected) {
            expect_identical(
                object = unname(convertUCSCBuildToEnsembl(object)),
                expected = expected
            )
        },
        SIMPLIFY = FALSE
    )
})

test_that("Expected failure", {
    expect_error(
        object = convertUCSCBuildToEnsembl("XXX"),
        regexp = "Failed to match UCSC"
    )
})
