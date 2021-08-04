context("downloadGenome")

test_that("downloadEnsemblGenome", {
    info <- downloadEnsemblGenome(
        organism = "Homo sapiens",
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

test_that("downloadGencodeGenome", {
    info <- downloadGencodeGenome(
        organism = "Homo sapiens",
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
