test_that("ncbiGeneId (IntegerList) serializes as comma-separated integers", {
    ## Verify that CompressedIntegerList columns export via toString() to
    ## comma-separated character when writing to CSV — matching the behaviour
    ## of NCBI Gene ID mappings on EnsemblGenes objects.
    object <- GRanges(
        seqnames = "1",
        ranges = IRanges::IRanges(start = c(1L, 100L), width = 10L),
        strand = "+"
    )
    mcols(object)[["geneId"]] <- c("ENSG0001", "ENSG0002")
    mcols(object)[["ncbiGeneId"]] <- IntegerList(
        list(4780L, c(1234L, 5678L))
    )
    names(object) <- mcols(object)[["geneId"]]
    con <- file.path(tempdir2(), "genes_ncbi.csv")
    on.exit(unlink2(con))
    export(object = object, con = con)
    expect_true(file.exists(con))
    df <- import(con = con)
    expect_true(isSubset("ncbiGeneId", colnames(df)))
    ## Single NCBI Gene ID should be character "4780".
    expect_identical(df[["ncbiGeneId"]][[1L]], "4780")
    ## Multiple NCBI Gene IDs should be comma-separated "1234, 5678".
    expect_identical(df[["ncbiGeneId"]][[2L]], "1234, 5678")
})


test_that("TxToGene", {
    object <- DataFrame(
        txId = c(
            "tx0001",
            "tx0002",
            "tx0003",
            "tx0004"
        ),
        geneId = c(
            "gene0001",
            "gene0001",
            "gene0002",
            "gene0002"
        )
    )
    object <- TxToGene(object)
    expect_s4_class(object, "TxToGene")
    con <- file.path(tempdir2(), "tx2gene.csv")
    expect_false(file.exists(con))
    export(object = object, con = con)
    expect_true(file.exists(con))
    expect_identical(
        object = import(
            con = con,
            format = "lines",
            nMax = 4L
        ),
        expected = c(
            "\"tx0001\",\"gene0001\"",
            "\"tx0002\",\"gene0001\"",
            "\"tx0003\",\"gene0002\"",
            "\"tx0004\",\"gene0002\""
        )
    )
    unlink2(con)
})
