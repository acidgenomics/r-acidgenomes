context("matchesGene2Symbol")

gene2symbol <- Gene2Symbol(
    object = DataFrame(
        geneID = paste0("gene", seq_len(2L)),
        geneName = paste0("symbol", seq_len(2L))
    )
)

test_that("matchesGene2Symbol", {
    genes <- gene2symbol[["geneID"]]
    expect_true(is.character(genes))
    x <- DataFrame(
        "sample1" = c(1L, 2L),
        "sample2" = c(3L, 4L),
        row.names = genes
    )
    expect_true(
        matchesGene2Symbol(
            x = x,
            genes = genes,
            gene2symbol = gene2symbol
        )
    )
    expect_false(
        matchesGene2Symbol(
            x = datasets::mtcars,
            genes = genes,
            gene2symbol = gene2symbol
        )
    )
})
