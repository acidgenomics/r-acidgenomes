context("downloadGenome")

test_that("downloadEnsemblGenome", {
    info <- downloadEnsemblGenome(
        organism = "Homo sapiens",
        genomeBuild = "GRCh38",
        release = 104L,
        cache = TRUE
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
    expect_identical(nrow(tx2gene), 257575L)
    expect_identical(
        object = as.data.frame(tx2gene)[1L, , drop = TRUE],
        expected = list(
            "txId" = "ENST00000000233.10",
            "geneId" = "ENSG00000004059.11"
        )
    )
    tx2gene <- import(
        file = file.path(outputDir, "tx2gene.csv.gz"),
        colnames = FALSE
    )
    expect_identical(nrow(tx2gene), 257575L)
    expect_identical(
        object = as.data.frame(tx2gene)[1L, , drop = TRUE],
        expected = list(
            "V1" = "ENST00000000233.10",
            "V2" = "ENSG00000004059.11"
        )
    )
    if (dir.exists(outputDir)) {
        unlink(outputDir, recursive = TRUE)
    }
})

test_that("downloadGencodeGenome", {
    info <- downloadGencodeGenome(
        organism = "Homo sapiens",
        genomeBuild = "GRCh38",
        release = 38L,
        cache = TRUE
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
    expect_identical(nrow(tx2gene), 237012L)
    expect_identical(
        object = as.data.frame(tx2gene)[1L, , drop = TRUE],
        expected = list(
            "txId" = "ENST00000000233.10",
            "geneId" = "ENSG00000004059.11"
        )
    )
    tx2gene <- import(
        file = file.path(outputDir, "tx2gene.csv.gz"),
        colnames = FALSE
    )
    expect_identical(nrow(tx2gene), 237012L)
    expect_identical(
        object = as.data.frame(tx2gene)[1L, , drop = TRUE],
        expected = list(
            "V1" = "ENST00000000233.10",
            "V2" = "ENSG00000004059.11"
        )
    )
    if (dir.exists(outputDir)) {
        unlink(outputDir, recursive = TRUE)
    }
})

test_that("downloadRefSeqGenome", {
    info <- downloadRefSeqGenome(
        organism = "Homo sapiens",
        taxonomicGroup = "vertebrate_mammalian",
        cache = TRUE
    )
    outputDir <- info[["args"]][["outputDir"]]
    expect_true(dir.exists(outputDir))
    expect_true(all(file.exists(file.path(
        outputDir,
        c(
            "annotation.gff3.gz",
            "annotation.gtf.gz",
            "genome.fa.gz",
            "metadata.rds",
            "transcriptome.fa.gz",
            "tx2gene.csv.gz"
        )
    ))))
    if (dir.exists(outputDir)) {
        unlink(outputDir, recursive = TRUE)
    }
})

test_that("downloadUCSCGenome", {
    info <- downloadUCSCGenome(
        organism = "Homo sapiens",
        cache = TRUE
    )
    outputDir <- info[["args"]][["outputDir"]]
    expect_true(dir.exists(outputDir))
    expect_true(all(file.exists(file.path(
        outputDir,
        c(
            "annotation.gtf.gz",
            "genome.fa.gz",
            "metadata.rds",
            "transcriptome.fa.gz",
            "tx2gene.csv.gz"
        )
    ))))
    if (dir.exists(outputDir)) {
        unlink(outputDir, recursive = TRUE)
    }
})
