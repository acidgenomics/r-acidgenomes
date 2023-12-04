test_that("Homo sapiens GRCh38", {
    testdir <- tempdir2()
    info <- downloadRefseqGenome(
        organism = "Homo sapiens",
        taxonomicGroup = "vertebrate_mammalian",
        genomeBuild = "GCF_000001405.40_GRCh38.p14",
        outputDir = testdir,
        cache = TRUE
    )
    outputDir <- info[["args"]][["outputDir"]]
    expect_true(dir.exists(outputDir))
    expect_match(
        object = basename(outputDir),
        regexp = "homo-sapiens-gcf-000001405-40-grch38-p14-refseq"
    )
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
    expect_s4_class(genes, "RefseqGenes")
    expect_identical(
        object = head(sort(names(genes)), n = 3L),
        expected = c("A1BG", "A1BG-AS1", "A1CF")
    )
    transcripts <- import(file.path(outputDir, "transcripts.rds"))
    expect_s4_class(transcripts, "RefseqTranscripts")
    expect_identical(
        object = head(sort(names(transcripts)), n = 3L),
        expected = c("NM_000014.6", "NM_000015.3", "NM_000016.6")
    )
    t2g <- import(file.path(outputDir, "tx2gene.rds"))
    expect_s4_class(t2g, "TxToGene")
    t2g <- as.data.frame(t2g)
    aatfExpected <- data.frame(
        "txId" = c(
            "NM_001411094.1",
            "NM_012138.4",
            "XM_047435748.1",
            "XM_054329280.1"
        ),
        "geneId" = rep("AATF", 4L)
    )
    aatfCurrent <- t2g[t2g[, 2L] == "AATF", ]
    rownames(aatfCurrent) <- NULL
    expect_identical(aatfCurrent, aatfExpected)
    t2g <- import(
        con = file.path(outputDir, "tx2gene.csv.gz"),
        colnames = c("txId", "geneId")
    )
    aatfCurrent <- t2g[t2g[, 2L] == "AATF", ]
    rownames(aatfCurrent) <- NULL
    expect_identical(aatfCurrent, aatfExpected)
    unlink2(testdir)
})

test_that("Homo sapiens GRCh37", {
    testdir <- tempdir2()
    info <- downloadRefseqGenome(
        organism = "Homo sapiens",
        taxonomicGroup = "vertebrate_mammalian",
        genomeBuild = "GCF_000001405.25_GRCh37.p13",
        outputDir = testdir,
        cache = TRUE
    )
    outputDir <- info[["args"]][["outputDir"]]
    expect_true(dir.exists(outputDir))
    expect_match(
        object = basename(outputDir),
        regexp = "homo-sapiens-gcf-000001405-25-grch37-p13-refseq"
    )
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
    expect_s4_class(genes, "RefseqGenes")
    expect_identical(
        object = head(sort(names(genes)), n = 3L),
        expected = c("A1BG", "A1BG-AS1", "A1CF")
    )
    transcripts <- import(file.path(outputDir, "transcripts.rds"))
    expect_s4_class(transcripts, "RefseqTranscripts")
    expect_identical(
        object = head(sort(names(transcripts)), n = 3L),
        expected = c("NM_000014.6", "NM_000015.3", "NM_000016.6")
    )
    t2g <- import(file.path(outputDir, "tx2gene.rds"))
    expect_s4_class(t2g, "TxToGene")
    t2g <- as.data.frame(t2g)
    aatfExpected <- data.frame(
        "txId" = c("NM_012138.4"),
        "geneId" = "AATF"
    )
    aatfCurrent <- t2g[t2g[, 2L] == "AATF", ]
    rownames(aatfCurrent) <- NULL
    expect_identical(aatfCurrent, aatfExpected)
    t2g <- import(
        con = file.path(outputDir, "tx2gene.csv.gz"),
        colnames = c("txId", "geneId")
    )
    aatfCurrent <- t2g[t2g[, 2L] == "AATF", ]
    rownames(aatfCurrent) <- NULL
    expect_identical(aatfCurrent, aatfExpected)
    unlink2(testdir)
})

test_that("Mus musculus", {
    testdir <- tempdir2()
    info <- downloadRefseqGenome(
        organism = "Mus musculus",
        taxonomicGroup = "vertebrate_mammalian",
        genomeBuild = "GCF_000001635.27_GRCm39",
        outputDir = testdir,
        cache = TRUE
    )
    outputDir <- info[["args"]][["outputDir"]]
    expect_true(dir.exists(outputDir))
    expect_match(
        object = basename(outputDir),
        regexp = "mus-musculus-gcf-000001635-27-grcm39-refseq"
    )
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
    expect_s4_class(genes, "RefseqGenes")
    expect_identical(
        object = head(sort(names(genes)), n = 3L),
        expected = c("0610005C13Rik", "0610006L08Rik", "0610009B22Rik")
    )
    transcripts <- import(file.path(outputDir, "transcripts.rds"))
    expect_s4_class(transcripts, "RefseqTranscripts")
    expect_identical(
        object = head(sort(names(transcripts)), n = 3L),
        expected = c("NM_001001130.3", "NM_001001144.3", "NM_001001152.2")
    )
    t2g <- import(file.path(outputDir, "tx2gene.rds"))
    expect_s4_class(t2g, "TxToGene")
    t2g <- as.data.frame(t2g)
    nfe2l2Expected <- data.frame(
        "txId" = c("NM_001399226.1", "NM_010902.5"),
        "geneId" = rep("Nfe2l2", 2L)
    )
    nfe2l2Current <- t2g[t2g[, 2L] == "Nfe2l2", ]
    rownames(nfe2l2Current) <- NULL
    expect_identical(nfe2l2Current, nfe2l2Expected)
    t2g <- import(
        con = file.path(outputDir, "tx2gene.csv.gz"),
        colnames = c("txId", "geneId")
    )
    nfe2l2Current <- t2g[t2g[, 2L] == "Nfe2l2", ]
    rownames(nfe2l2Current) <- NULL
    expect_identical(nfe2l2Current, nfe2l2Expected)
    unlink2(testdir)
})
