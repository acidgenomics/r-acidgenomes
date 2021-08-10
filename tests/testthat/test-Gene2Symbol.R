context("Gene2Symbol")

test_that("Gene2Symbol", {
    formats <- eval(formals(`Gene2Symbol,GRanges`)[["format"]])
    for (format in formats) {
        object <- gr
        object <- Gene2Symbol(object, format = format)
        expect_s4_class(object, "Gene2Symbol")
        expect_identical(colnames(object), c("geneId", "geneName"))
    }
})

test_that("makeUnique mode", {
    object <- Gene2Symbol(
        object = DataFrame(
            "geneId" = c(
                "gene1",
                "gene2",
                "gene3",
                "gene4"
            ),
            "geneName" = c(
                "symbol1",
                "symbol1",
                "symbol2",
                NA_character_
            )
        ),
        format = "makeUnique"
    )
    expected <- DataFrame(
        "geneId" = c(
            "gene1",
            "gene2",
            "gene3",
            "gene4"
        ),
        "geneName" = c(
            "symbol1",
            "symbol1.1",
            "symbol2",
            "unannotated"
        )
    )
    metadata(expected) <- list(
        "format" = "makeUnique",
        "dupes" = "symbol1"
    )
    expected <- new(Class = "Gene2Symbol", expected)
    expect_identical(object, expected)
    ## Don't allow "geneId" column to contain NA values.
    expect_error(
        Gene2Symbol(
            object = DataFrame(
                "geneId" = c(
                    "gene1",
                    NA_character_
                ),
                "geneName" = c(
                    "symbol1",
                    "symbol2"
                )
            ),
            format = "makeUnique"
        )
    )
    ## Ensure gene identifiers return sorted.
    object <- Gene2Symbol(
        object = DataFrame(
            "geneId" = c(
                "gene2",
                "gene1"
            ),
            "geneName" = c(
                "symbol2",
                "symbol1"
            ),
            row.names = c("B", "A")
        ),
        format = "makeUnique"
    )
    expected <- DataFrame(
        "geneId" = c("gene1", "gene2"),
        "geneName" = c("symbol1", "symbol2"),
        row.names = c("A", "B")
    )
    metadata(expected) <- list("format" = "makeUnique")
    expected <- new(Class = "Gene2Symbol", expected)
    expect_identical(object, expected)
})

test_that("1:1 mode", {
})

test_that("unmodified mode", {
})

test_that("summary", {
    object <- gr
    x <- Gene2Symbol(object)
    output <- capture.output(summary(x))
    expect_identical(
        head(output, n = 1L),
        paste("genes:", length(object))
    )
})
