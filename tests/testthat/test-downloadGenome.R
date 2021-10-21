context("downloadGenome")

testdir <- file.path(tempdir(), "genome")

## FIXME This may not be saving tx2gene to disk correctly, resulting in the
## number of transcript-to-gene rows being off by 1....argh.
## FIXME This likely is due to the BiocIO export migration...

test_that("downloadEnsemblGenome", {
    unlink(testdir, recursive = TRUE)
    info <- downloadEnsemblGenome(
        organism = "Homo sapiens",
        genomeBuild = "GRCh38",
        release = 104L,
        outputDir = testdir,
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


    ## FIXME Check the number of genes in GRanges.
    ## FIXME Check the number of transcripts in GRanges.

    tx2gene <- import(file.path(outputDir, "tx2gene.rds"))
    expect_identical(nrow(tx2gene), 257575L)
    expect_identical(
        object = as.data.frame(tx2gene)[1L, , drop = TRUE],
        expected = list(
            "txId" = "ENST00000000233.10",
            "geneId" = "ENSG00000004059.11"
        )
    )

    ## FIXME This method isn't inheriting in package correctly...
    ## FIXME This BiocIO change is problematic with file, need to rework...
    ## FIXME The export method is including colnames for tx2gene.csv....
    ## we don't want this behavior...
    ## FIXME Double check that tx2gene.csv DOESN'T contain column names.

    tx2gene <- import(
        file = file.path(outputDir, "tx2gene.csv.gz"),
        colnames = c("txId", "geneId")
    )

    expect_identical(nrow(tx2gene), 257575L)
    expect_identical(
        object = as.data.frame(tx2gene)[1L, , drop = TRUE],
        expected = list(
            "txId" = "ENST00000000233.10",
            "geneId" = "ENSG00000004059.11"
        )
    )
    unlink(testdir, recursive = TRUE)
})

test_that("downloadGencodeGenome", {
    unlink(testdir, recursive = TRUE)
    info <- downloadGencodeGenome(
        organism = "Homo sapiens",
        genomeBuild = "GRCh38",
        release = 38L,
        outputDir = testdir,
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
        colnames = c("txId", "geneId")
    )
    expect_identical(nrow(tx2gene), 237012L)
    expect_identical(
        object = as.data.frame(tx2gene)[1L, , drop = TRUE],
        expected = list(
            "txId" = "ENST00000000233.10",
            "geneId" = "ENSG00000004059.11"
        )
    )
    unlink(testdir, recursive = TRUE)
})

test_that("downloadRefSeqGenome", {
    unlink(testdir, recursive = TRUE)
    info <- downloadRefSeqGenome(
        organism = "Homo sapiens",
        taxonomicGroup = "vertebrate_mammalian",
        genomeBuild = "GCF_000001405.39_GRCh38.p13",
        outputDir = testdir,
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
    unlink(testdir, recursive = TRUE)
})

test_that("downloadUCSCGenome", {
    unlink(testdir, recursive = TRUE)
    info <- downloadUCSCGenome(
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
    unlink(testdir, recursive = TRUE)
})
