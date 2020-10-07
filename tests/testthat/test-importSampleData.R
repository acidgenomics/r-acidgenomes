context("importSampleData : Demultiplexed samples")

file <- file.path("cache", "bcbio-metadata-demultiplexed.csv")

test_that("DataFrame return", {
    expect_identical(
        object = importSampleData(file, pipeline = "bcbio"),
        expected = DataFrame(
            sampleName = factor(paste0("sample", seq_len(4L))),
            fileName = factor(paste0("sample", seq_len(4L), "_R1.fastq.gz")),
            description = factor(paste0("sample", seq_len(4L))),
            genotype = factor(rep(c("wildtype", "knockout"), times = 2L)),
            row.names = paste0("sample", seq_len(4L))
        )
    )
})

test_that("Lane-split technical replicate support", {
    object <- importSampleData(file, lanes = 4L, pipeline = "bcbio")
    expect_true("lane" %in% colnames(object))
    expect_identical(
        object = rownames(object)[seq_len(8L)],
        expected = c(
            paste0("sample1_L00", seq_len(4L)),
            paste0("sample2_L00", seq_len(4L))
        )
    )
})

test_that("Required column check failure", {
    file <- file.path(
        "cache",
        "bcbio-metadata-demultiplexed-invalid-missing-columns.csv"
    )
    expect_error(
        object = importSampleData(file, pipeline = "bcbio"),
        regexp = "description"
    )
})

test_that("Duplicated description", {
    file <- file.path(
        "cache",
        "bcbio-metadata-demultiplexed-invalid-duplicated.csv"
    )
    expect_error(
        object = importSampleData(file, pipeline = "bcbio"),
        regexp = "Sample data input file is malformed."
    )
})



context("importSampleData : Multiplexed samples")

file <- file.path("cache", "bcbio-metadata-multiplexed-indrops.csv")

test_that("DataFrame return", {
    object <- importSampleData(file, pipeline = "bcbio")
    expected <- DataFrame(
        sampleName = factor(c(
            "sample2_1",
            "sample1_1",
            "sample3_1",
            "sample4_1",
            "sample2_2",
            "sample1_2",
            "sample3_2",
            "sample4_2"
        )),
        fileName = factor(c(
            rep("indrops1_R1.fastq.gz", times = 4L),
            rep("indrops2_R1.fastq.gz", times = 4L)
        )),
        ## Valid rownames (sampleIDs) are generated from this column.
        ## Note that we're sorting the sample metadata by this column.
        description = factor(c(
            "indrops1-AGAGGATA",
            "indrops1-ATAGAGAG",
            "indrops1-CTCCTTAC",
            "indrops1-TATGCAGT",
            "indrops2-AGAGGATA",
            "indrops2-ATAGAGAG",
            "indrops2-CTCCTTAC",
            "indrops2-TATGCAGT"
        )),
        index = factor(rep(c(2L, 1L, 3L, 4L), times = 2L)),
        sequence = factor(
            rep(c(
                "TATCCTCT",
                "CTCTCTAT",
                "GTAAGGAG",
                "ACTGCATA"
            ), times = 2L)
        ),
        aggregate = factor(
            paste0("sample", rep(c(2L, 1L, 3L, 4L), times = 2L)),
            levels = paste0("sample", seq_len(4L))
        ),
        genotype = factor(
            rep(
                x = c(
                    "knockout",
                    "wildtype",
                    "wildtype",
                    "knockout"
                ),
                times = 2L
            ),
            ## Note that the order should be alphabetical here.
            levels = c("knockout", "wildtype")
        ),
        revcomp = factor(
            rep(c(
                "AGAGGATA",
                "ATAGAGAG",
                "CTCCTTAC",
                "TATGCAGT"
            ), times = 2L)
        ),
        row.names = c(
            "indrops1_AGAGGATA",
            "indrops1_ATAGAGAG",
            "indrops1_CTCCTTAC",
            "indrops1_TATGCAGT",
            "indrops2_AGAGGATA",
            "indrops2_ATAGAGAG",
            "indrops2_CTCCTTAC",
            "indrops2_TATGCAGT"
        )
    )
    expect_identical(object, expected)
})

test_that("Lane-split technical replicate support", {
    object <- importSampleData(file, lanes = 4L, pipeline = "bcbio")
    expect_identical(
        object = rownames(object),
        expected = c(
            "indrops1_L001_AGAGGATA",
            "indrops1_L001_ATAGAGAG",
            "indrops1_L001_CTCCTTAC",
            "indrops1_L001_TATGCAGT",
            "indrops1_L002_AGAGGATA",
            "indrops1_L002_ATAGAGAG",
            "indrops1_L002_CTCCTTAC",
            "indrops1_L002_TATGCAGT",
            "indrops1_L003_AGAGGATA",
            "indrops1_L003_ATAGAGAG",
            "indrops1_L003_CTCCTTAC",
            "indrops1_L003_TATGCAGT",
            "indrops1_L004_AGAGGATA",
            "indrops1_L004_ATAGAGAG",
            "indrops1_L004_CTCCTTAC",
            "indrops1_L004_TATGCAGT",
            "indrops2_L001_AGAGGATA",
            "indrops2_L001_ATAGAGAG",
            "indrops2_L001_CTCCTTAC",
            "indrops2_L001_TATGCAGT",
            "indrops2_L002_AGAGGATA",
            "indrops2_L002_ATAGAGAG",
            "indrops2_L002_CTCCTTAC",
            "indrops2_L002_TATGCAGT",
            "indrops2_L003_AGAGGATA",
            "indrops2_L003_ATAGAGAG",
            "indrops2_L003_CTCCTTAC",
            "indrops2_L003_TATGCAGT",
            "indrops2_L004_AGAGGATA",
            "indrops2_L004_ATAGAGAG",
            "indrops2_L004_CTCCTTAC",
            "indrops2_L004_TATGCAGT"
        )
    )
})

test_that("Required column check failure.", {
    file <- file.path(
        "cache",
        "bcbio-metadata-multiplexed-invalid-missing-columns.csv"
    )
    expect_error(
        object = importSampleData(file, pipeline = "bcbio"),
        regexp = "Sample data input file is malformed."
    )
})

test_that("Duplicate rows in 'sampleName' column", {
    file <- file.path(
        "cache",
        "bcbio-metadata-multiplexed-invalid-duplicated.csv"
    )
    expect_error(
        object = importSampleData(file, pipeline = "bcbio"),
        regexp = "sampleName"
    )
})

## Recommend using `fileName` instead.
test_that("bcbio 'samplename' column", {
    file <- file.path(
        "cache",
        "bcbio-metadata-demultiplexed-invalid-legacy-samplename.csv"
    )
    out <- importSampleData(file, pipeline = "bcbio")
    expect_identical(
        object = colnames(out),
        expected = c("sampleName", "fileName", "description")
    )
})

test_that("'sampleID' column defined by user", {
    file <- file.path(
        "cache",
        "bcbio-metadata-demultiplexed-invalid-sample-id.csv"
    )
    expect_error(
        object = importSampleData(file, pipeline = "bcbio"),
        regexp = "sampleID"
    )
    expect_s4_class(
        object = importSampleData(file, pipeline = "none"),
        class = "DataFrame"
    )
})

test_that("Missing file", {
    expect_error(
        object = importSampleData("XXX.csv"),
        regexp = "isAFile"
    )
})



context("importSampleData : Malformed input")

test_that("Metadata blacklist", {
    file <- file.path(
        "cache",
        "bcbio-metadata-invalid-column-name.csv"
    )
    expect_error(
        object = importSampleData(file, pipeline = "bcbio"),
        regexp = "sampleNames"
    )
})

test_that("Invalid description", {
    file <- file.path(
        "cache",
        "bcbio-metadata-invalid-description.csv"
    )
    expect_error(
        object = importSampleData(file, pipeline = "bcbio"),
        regexp = "Sample data input file is malformed."
    )
})
