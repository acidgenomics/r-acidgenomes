test_that("downloadEnsemblGenome", {
    testdir <- tempdir2()
    info <- downloadEnsemblGenome(
        organism = "Homo sapiens",
        genomeBuild = "GRCh38",
        release = 105L,
        outputDir = testdir,
        cache = FALSE
    )
    outputDir <- info[["args"]][["outputDir"]]
    expect_true(dir.exists(outputDir))
    expect_true(all(file.exists(file.path(
        outputDir,
        c(
            "annotation.gff3.gz",
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
    expect_s4_class(genes, "EnsemblGenes")
    expect_length(genes, 61541L)
    expect_identical(
        object = head(sort(names(genes)), n = 3L),
        expected = c(
            "ENSG00000000003.15",
            "ENSG00000000005.6",
            "ENSG00000000419.14"
        )
    )
    transcripts <- import(file.path(outputDir, "transcripts.rds"))
    expect_length(transcripts, 244825L)
    expect_identical(
        object = head(sort(names(transcripts)), n = 3L),
        expected = c(
            "ENST00000000233.10",
            "ENST00000000412.8",
            "ENST00000000442.11"
        )
    )
    tx2gene <- import(file.path(outputDir, "tx2gene.rds"))
    expect_identical(nrow(tx2gene), 265480L)
    expect_identical(
        object = as.data.frame(tx2gene)[1L, , drop = TRUE],
        expected = list(
            "txId" = "ENST00000000233.10",
            "geneId" = "ENSG00000004059.11"
        )
    )
    tx2gene <- import(
        file = file.path(outputDir, "tx2gene.csv.gz"),
        colnames = c("txId", "geneId")
    )
    expect_identical(nrow(tx2gene), 265480L)
    expect_identical(
        object = as.data.frame(tx2gene)[1L, , drop = TRUE],
        expected = list(
            "txId" = "ENST00000000233.10",
            "geneId" = "ENSG00000004059.11"
        )
    )
    unlink2(testdir)
})

test_that("downloadGencodeGenome", {
    testdir <- tempdir2()
    info <- downloadGencodeGenome(
        organism = "Homo sapiens",
        genomeBuild = "GRCh38",
        release = 39L,
        outputDir = testdir,
        cache = FALSE
    )
    outputDir <- info[["args"]][["outputDir"]]
    expect_true(dir.exists(outputDir))
    expect_true(all(file.exists(file.path(
        outputDir,
        c(
            "annotation.gff3.gz",
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
    tx2gene <- import(file.path(outputDir, "tx2gene.rds"))
    expect_identical(nrow(tx2gene), 244939L)
    expect_identical(
        object = as.data.frame(tx2gene)[1L, , drop = TRUE],
        expected = list(
            "txId" = "ENST00000000233.10",
            "geneId" = "ENSG00000004059.11"
        )
    )
    tx2gene <- import(
        file = file.path(outputDir, "tx2gene.csv.gz"),
        colnames = c("txId", "geneId")
    )
    expect_identical(nrow(tx2gene), 244939L)
    expect_identical(
        object = as.data.frame(tx2gene)[1L, , drop = TRUE],
        expected = list(
            "txId" = "ENST00000000233.10",
            "geneId" = "ENSG00000004059.11"
        )
    )
    unlink2(testdir)
})

test_that("downloadRefSeqGenome", {
    testdir <- tempdir2()
    info <- downloadRefSeqGenome(
        organism = "Homo sapiens",
        taxonomicGroup = "vertebrate_mammalian",
        genomeBuild = "GCF_000001405.39_GRCh38.p13",
        outputDir = testdir,
        cache = FALSE
    )
    outputDir <- info[["args"]][["outputDir"]]
    expect_true(dir.exists(outputDir))
    expect_true(all(file.exists(file.path(
        outputDir,
        c(
            "annotation.gff3.gz",
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
    aatfExpected <- data.frame(
        "txId" = c(
            "NM_012138.4",
            "XM_011524611.2",
            "XR_934439.3"
        ),
        "geneId" = rep("AATF", 3L)
    )
    tx2gene <- import(file.path(outputDir, "tx2gene.rds"))
    tx2gene <- as.data.frame(tx2gene)
    aatfCurrent <- tx2gene[tx2gene[, 2L] == "AATF", ]
    rownames(aatfCurrent) <- NULL
    expect_identical(aatfCurrent, aatfExpected)
    tx2gene <- import(
        file = file.path(outputDir, "tx2gene.csv.gz"),
        colnames = c("txId", "geneId")
    )
    aatfCurrent <- tx2gene[tx2gene[, 2L] == "AATF", ]
    rownames(aatfCurrent) <- NULL
    expect_identical(aatfCurrent, aatfExpected)
    unlink2(testdir)
})

test_that("downloadUCSCGenome", {
    testdir <- tempdir2()
    info <- downloadUCSCGenome(
        organism = "Homo sapiens",
        genomeBuild = "hg38",
        outputDir = testdir,
        cache = FALSE
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
    tx2gene <- import(file.path(outputDir, "tx2gene.rds"))
    expect_identical(
        object = as.data.frame(tx2gene)[1L, , drop = TRUE],
        expected = list(
            "txId" = "ENST00000000233",
            "geneId" = "ENSG00000004059"
        )
    )
    tx2gene <- import(
        file = file.path(outputDir, "tx2gene.csv.gz"),
        colnames = c("txId", "geneId")
    )
    expect_identical(
        object = as.data.frame(tx2gene)[1L, , drop = TRUE],
        expected = list(
            "txId" = "ENST00000000233",
            "geneId" = "ENSG00000004059"
        )
    )
    unlink2(testdir)
})
