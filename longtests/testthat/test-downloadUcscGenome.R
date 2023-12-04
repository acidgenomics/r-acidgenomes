## FIXME Add coverage for hg19.

test_that("downloadUcscGenome", {
    testdir <- tempdir2()
    info <- downloadUcscGenome(
        organism = "Homo sapiens",
        genomeBuild = "hg38",
        outputDir = testdir,
        cache = TRUE
    )
    outputDir <- info[["args"]][["outputDir"]]
    expect_true(dir.exists(outputDir))
    expect_true(all(file.exists(file.path(
        outputDir,
        c(
            "annotation.gtf.gz",
            "genes.rds",
            "genome.fa.gz",
            "metadata.rds",
            "transcriptome.fa.gz",
            "transcripts.rds",
            "tx2gene.csv.gz",
            "tx2gene.rds"
        )
    ))))
    genes <- import(file.path(outputDir, "genes.rds"))
    expect_s4_class(genes, "UcscGenes")
    expect_identical(
        object = head(sort(names(genes)), n = 3L),
        expected = c("A1BG", "A1BG-AS1", "A1CF")
    )
    transcripts <- import(file.path(outputDir, "transcripts.rds"))
    expect_s4_class(transcripts, "UcscTranscripts")
    expect_identical(
        object = head(sort(names(transcripts)), n = 3L),
        expected = c("ABCA3P1_2", "ABCB10P1_2", "ABCB10P1_3")
    )
    t2g <- import(file.path(outputDir, "tx2gene.rds"))
    expect_s4_class(t2g, "TxToGene")
    expect_identical(
        object = as.data.frame(t2g)[1L, , drop = TRUE],
        expected = list(
            "txId" = "ABCA3P1_2",
            "geneId" = "ABCA3P1"
        )
    )
    t2g <- import(
        con = file.path(outputDir, "tx2gene.csv.gz"),
        colnames = c("txId", "geneId")
    )
    expect_identical(
        object = as.data.frame(t2g)[1L, , drop = TRUE],
        expected = list(
            "txId" = "ABCA3P1_2",
            "geneId" = "ABCA3P1"
        )
    )
    unlink2(testdir)
})
