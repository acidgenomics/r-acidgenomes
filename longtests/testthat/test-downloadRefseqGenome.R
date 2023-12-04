## FIXME Add coverage for GRCh37.
## FIXME Add coverage for mouse genome.

## NOTE Add coverage for GRCh37/hg19 in a future update:
## https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/
## Homo_sapiens/all_assembly_versions/GCF_000001405.25_GRCh37.p13/

test_that("downloadRefseqGenome", {
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
