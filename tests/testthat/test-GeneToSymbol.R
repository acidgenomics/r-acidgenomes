test_that("GRanges format support", {
    formats <- eval(formals(`GeneToSymbol,DFrame`)[["format"]])
    for (format in formats) {
        object <- gr
        object <- GeneToSymbol(object, format = format)
        expect_s4_class(object, "GeneToSymbol")
    }
})

test_that("format : makeUnique", {
    object <- GeneToSymbol(
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
            "gene3"
        ),
        "geneName" = c(
            "symbol1",
            "symbol1.1",
            "symbol2"
        )
    )
    metadata(expected) <- list(
        "date" = Sys.Date(),
        "dropped" = "gene4",
        "dupes" = "symbol1",
        "format" = "makeUnique",
        "packageVersion" = .pkgVersion
    )
    expected <- new(Class = "GeneToSymbol", expected)
    expect_identical(object, expected)
    expect_error(
        object = GeneToSymbol(
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
        ),
        regexp = "anyNA"
    )
    ## Ensure gene identifiers return sorted.
    object <- GeneToSymbol(
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
    metadata(expected) <- list(
        "date" = Sys.Date(),
        "format" = "makeUnique",
        "packageVersion" = .pkgVersion
    )
    expected <- new(Class = "GeneToSymbol", expected)
    expect_identical(object, expected)
})

test_that("format : 1:1", {
    object <- GeneToSymbol(
        object = DataFrame(
            "geneId" = c(
                "gene2",
                "gene1",
                "gene4",
                "gene3",
                "gene5",
                "gene6"
            ),
            "geneName" = c(
                "symbol1",
                "symbol1",
                "symbol2",
                "symbol2",
                NA_character_,
                "symbol3"
            ),
            row.names = LETTERS[seq_len(6L)]
        ),
        format = "1:1"
    )
    expected <- DataFrame(
        "geneId" = c(
            "gene1",
            "gene3",
            "gene6"
        ),
        "geneName" = c(
            "symbol1",
            "symbol2",
            "symbol3"
        ),
        row.names = c("B", "D", "F")
    )
    metadata(expected) <- list(
        "date" = Sys.Date(),
        "dropped" = "gene5",
        "dupes" = c("symbol1", "symbol2"),
        "format" = "1:1",
        "packageVersion" = .pkgVersion
    )
    expected <- new(Class = "GeneToSymbol", expected)
    expect_identical(object, expected)
})

test_that("format : unmodified", {
    object <- GeneToSymbol(
        object = DataFrame(
            "geneId" = c(
                "gene1",
                "gene2",
                "gene3",
                "gene4",
                "gene5"
            ),
            "geneName" = c(
                "symbol1",
                "symbol1",
                "symbol2",
                NA_character_,
                "symbol3"
            ),
            row.names = LETTERS[seq_len(5L)]
        ),
        format = "unmodified"
    )
    expected <- DataFrame(
        "geneId" = c(
            "gene1",
            "gene2",
            "gene3",
            "gene5"
        ),
        "geneName" = c(
            "symbol1",
            "symbol1",
            "symbol2",
            "symbol3"
        ),
        row.names = c("A", "B", "C", "E")
    )
    metadata(expected) <- list(
        "date" = Sys.Date(),
        "dropped" = "gene4",
        "dupes" = "symbol1",
        "format" = "unmodified",
        "packageVersion" = .pkgVersion
    )
    expected <- new(Class = "GeneToSymbol", expected)
    expect_identical(object, expected)
})

test_that("summary method", {
    object <- gr
    x <- GeneToSymbol(object)
    output <- capture.output(summary(x))
    expect_identical(
        head(output, n = 1L),
        paste("genes:", length(object))
    )
})
