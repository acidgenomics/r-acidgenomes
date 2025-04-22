test_that("Homo sapiens hg38", {
    testdir <- tempdir2()
    info <- downloadUcscGenome(
        organism = "Homo sapiens",
        genomeBuild = "hg38",
        outputDir = testdir,
        cache = TRUE
    )
    outputDir <- info[["args"]][["outputDir"]]
    expect_true(dir.exists(outputDir))
    expect_match(
        object = basename(outputDir),
        regexp = "homo-sapiens-hg38-ucsc"
    )
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

test_that("Homo sapiens hg19", {
    testdir <- tempdir2()
    info <- downloadUcscGenome(
        organism = "Homo sapiens",
        genomeBuild = "hg19",
        outputDir = testdir,
        cache = TRUE
    )
    outputDir <- info[["args"]][["outputDir"]]
    expect_true(dir.exists(outputDir))
    expect_match(
        object = basename(outputDir),
        regexp = "homo-sapiens-hg19-ucsc"
    )
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
        expected = c("ABCB10P4_2", "ABCD1P4_2", "ABHD17AP1_2")
    )
    t2g <- import(file.path(outputDir, "tx2gene.rds"))
    expect_s4_class(t2g, "TxToGene")
    expect_identical(
        object = as.data.frame(t2g)[1L, , drop = TRUE],
        expected = list(
            "txId" = "ABCB10P4_2",
            "geneId" = "ABCB10P4"
        )
    )
    t2g <- import(
        con = file.path(outputDir, "tx2gene.csv.gz"),
        colnames = c("txId", "geneId")
    )
    expect_identical(
        object = as.data.frame(t2g)[1L, , drop = TRUE],
        expected = list(
            "txId" = "ABCB10P4_2",
            "geneId" = "ABCB10P4"
        )
    )
    unlink2(testdir)
})

test_that("Mus musculus mm39", {
    testdir <- tempdir2()
    ## Adding the collate step here to avoid different sorting of identifiers.
    with_collate(new = "C", code = {
        info <- downloadUcscGenome(
            organism = "Mus musculus",
            genomeBuild = "mm39",
            outputDir = testdir,
            cache = TRUE
        )
    })
    outputDir <- info[["args"]][["outputDir"]]
    expect_true(dir.exists(outputDir))
    expect_match(
        object = basename(outputDir),
        regexp = "mus-musculus-mm39-ucsc"
    )
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
        object = head(names(genes), n = 3L),
        expected = c("Gm26206", "Gm18956", "LOC118567655")
    )
    transcripts <- import(file.path(outputDir, "transcripts.rds"))
    expect_s4_class(transcripts, "UcscTranscripts")
    expect_identical(
        object = head(names(transcripts), n = 3L),
        expected = c("XR_004936710.1", "XR_004935518.1", "XR_004935517.1")
    )
    t2g <- import(file.path(outputDir, "tx2gene.rds"))
    expect_s4_class(t2g, "TxToGene")
    expect_identical(
        object = as.data.frame(t2g)[1L, , drop = TRUE],
        expected = list(
            "txId" = "NM_001001130.3",
            "geneId" = "Zfp85"
        )
    )
    t2g <- import(
        con = file.path(outputDir, "tx2gene.csv.gz"),
        colnames = c("txId", "geneId")
    )
    expect_identical(
        object = as.data.frame(t2g)[1L, , drop = TRUE],
        expected = list(
            "txId" = "NM_001001130.3",
            "geneId" = "Zfp85"
        )
    )
    unlink2(testdir)
})
