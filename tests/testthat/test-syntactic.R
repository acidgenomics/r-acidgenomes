data(syntactic, package = "AcidTest", envir = environment())
funs <- list(
    camelCase = camelCase,
    dottedCase = dottedCase,
    snakeCase = snakeCase,
    upperCamelCase = upperCamelCase
)



context("syntactic : atomic")

test_that("Unnamed", {
    for (f in funs) {
        object <- 1L
        expect_identical(f(object), object)
    }
})

test_that("Named", {
    mapply(
        f = funs,
        expected = c(
            "helloWorld",
            "hello.world",
            "hello_world",
            "HelloWorld"
        ),
        FUN = function(f, expected) {
            object <- c("hello.world" = 1L)
            expect_identical(names(f(object)), expected)
        },
        SIMPLIFY = FALSE
    )
})



context("syntactic : factor")

test_that("factor", {
    mapply(
        f = funs,
        levels = list(
            camelCase = c("group1", "group2"),
            dottedCase = c("group.1", "group.2"),
            snakeCase = c("group_1", "group_2"),
            upperCamelCase = c("Group1", "Group2")
        ),
        names = list(
            camelCase = c("sample1", "sample2", "sample3", "sample4"),
            dottedCase = c("sample.1", "sample.2", "sample.3", "sample.4"),
            snakeCase = c("sample_1", "sample_2", "sample_3", "sample_4"),
            upperCamelCase = c("Sample1", "Sample2", "Sample3", "Sample4")
        ),
        FUN = function(f, levels, names) {
            object <- syntactic[["factor"]]
            x <- f(object, names = TRUE)
            expect_identical(
                object = levels(x),
                expected = levels
            )
            expect_identical(
                object = names(x),
                expected = names
            )
            x <- f(object, names = FALSE)
            expect_identical(
                object = names(x),
                expected = names(object)
            )
        },
        SIMPLIFY = FALSE
    )
})



context("syntactic : list")

test_that("list", {
    mapply(
        f = funs,
        expected = list(
            camel = c("itemA", "itemB"),
            dotted = c("Item.A", "Item.B"),
            snake = c("item_a", "item_b"),
            upperCamel = c("ItemA", "ItemB")
        ),
        FUN = function(f, expected) {
            object <- syntactic[["list"]]
            expect_identical(
                object = names(f(object)),
                expected = expected
            )
        },
        SIMPLIFY = FALSE
    )
})



context("syntactic : matrix")

test_that("matrix", {
    mapply(
        f = funs,
        expected = list(
            camelCase = list(
                c("alabama", "alaska", "arizona"),
                c("murder", "assault", "urbanPop")
            ),
            dottedCase = list(
                c("Alabama", "Alaska", "Arizona"),
                c("Murder", "Assault", "Urban.Pop")
            ),
            snakeCase = list(
                c("alabama", "alaska", "arizona"),
                c("murder", "assault", "urban_pop")
            ),
            upperCamelCase = list(
                c("Alabama", "Alaska", "Arizona"),
                c("Murder", "Assault", "UrbanPop")
            )
        ),
        FUN = function(f, expected) {
            object <- syntactic[["matrix"]][seq_len(3L), seq_len(3L)]
            expect_identical(
                object = dimnames(f(object, rownames = TRUE, colnames = TRUE)),
                expected = expected
            )
        },
        SIMPLIFY = FALSE
    )
})



context("syntactic : Vector")

object <- as(df, "Vector")
names(object) <- toupper(names(object))
mcols(object) <- DataFrame(TEST = seq_len(length(object)))
metadata(object) <- list(TEST = "XXX")

test_that("SimpleList", {
    for (f in funs) {
        x <- f(
            object = object,
            names = TRUE,
            mcols = TRUE,
            metadata = TRUE
        )
        expect_s4_class(x, "SimpleList")
    }
})



context("syntactic : DataFrame")

mcols(df) <- DataFrame(TEST = seq_len(ncol(df)))
metadata(df) <- list(TEST = "XXX")

test_that("DataFrame", {
    for (f in funs) {
        x <- f(
            object = df,
            rownames = TRUE,
            colnames = TRUE,
            mcols = TRUE,
            metadata = TRUE
        )
        expect_s4_class(x, "DataFrame")
    }
})



context("syntactic : GRanges")

test_that("GRanges", {
    mapply(
        f = funs,
        expected = list(
            camelCase = c("geneID", "geneName"),
            dottedCase = c("gene.ID", "gene.Name"),
            snakeCase = c("gene_id", "gene_name"),
            upperCamelCase = c("GeneID", "GeneName")
        ),
        FUN = function(f, expected) {
            x <- f(
                object = gr,
                names = TRUE,
                mcols = TRUE
            )
            expect_s4_class(x, "GRanges")
            expect_identical(
                object = colnames(mcols(x)),
                expected = expected
            )
        },
        SIMPLIFY = FALSE
    )
})



context("syntactic : SummarizedExperiment")

test_that("SummarizedExperiment", {
    for (f in funs) {
        x <- f(
            object = rse,
            rownames = TRUE,
            colnames = TRUE,
            assayNames = TRUE,
            rowData = TRUE,
            colData = TRUE,
            metadata = TRUE
        )
        expect_s4_class(x, "SummarizedExperiment")
    }
})
