context("markdown")

test_that("SummarizedExperiment", {
    output <- capture.output(markdown(rse))
    expect_identical(
        object = output[[3L]],
        expected = "|         |condition |sampleName |interestingGroups |"
    )
})



context("markdownHeader")

test_that("default", {
    md <- markdownHeader("Header")
    expect_is(md, "knit_asis")
})

test_that("level", {
    md <- markdownHeader("Header", level = 4L)
    expect_identical(
        object = as.character(md),
        expected = "#### Header\n"  # nolint
    )
})

test_that("invalid level", {
    expect_error(
        object = markdownHeader("Header", level = 8L),
        expected = "Markdown supports header levels 1-7."
    )
})

test_that("tabset", {
    md <- markdownHeader("Header", tabset = TRUE)
    expect_identical(
        object = as.character(md),
        expected = "## Header {.tabset}\n"  # nolint
    )
})

test_that("asis", {
    expect_output(
        object = markdownHeader("Header", asis = TRUE),
        regexp = "\\n\\n## Header\\n"  # nolint
    )
})



context("markdownLink")

test_that("markdownLink", {
    expect_output(
        object = markdownLink(
            text = "R",
            url = "https://www.r-project.org",
            title = "The R Project for Statistical Computing"
        ),
        regexp = paste(
            "\\[R\\]\\(https://www.r-project.org\\):",
            "The R Project for Statistical Computing"
        )
    )
})



context("markdownList")

test_that("markdownList", {
    groceries <- c("milk", "eggs")
    md <- markdownList(groceries)
    expect_identical(
        class(md),
        "knit_asis"
    )
    md <- markdownList(groceries) %>%
        as.character()
    expect_identical(
        md,
        c("- milk\n", "- eggs\n")
    )
    md <- markdownList(groceries, ordered = TRUE) %>%
        as.character()
    expect_identical(
        md,
        c("1. milk\n", "2. eggs\n")
    )
    expect_output(
        markdownList(groceries, asis = TRUE),
        "\\n- milk\\n- eggs\\n"
    )
})



test_that("markdownPlots", {
    skip_if_not(hasInternet())
    skip_on_appveyor()
    plotlist <- readRDS(file = file.path("cache", "plotlist.rds"))
    output <- capture.output()
    expect_output(
        object = markdownPlots(list = plotlist),
        regexp = "## continuous"

    )
})



context("markdownTables")

test_that("Return unmodified in (interactive) R session", {
    x <- markdownTables(
        list = list(
            mpg = head(datasets::iris),
            mtcars = head(datasets::mtcars)
        ),
        captions = c("iris", "mtcars"),
        force = FALSE
    )
    expect_is(x, "list")
    expect_identical(names(x), c("mpg", "mtcars"))
})

test_that("Simulated knit session", {
    x <- markdownTables(
        list = list(
            head(datasets::iris),
            head(datasets::mtcars)
        ),
        captions = c("iris", "mtcars"),
        force = TRUE
    )
    expect_is(x, "knit_asis")

})
