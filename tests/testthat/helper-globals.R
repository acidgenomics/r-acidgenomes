data(
    DFrame,
    GRanges,
    RangedSummarizedExperiment,
    SingleCellExperiment,
    SummarizedExperiment_transcripts,
    matrix,
    matrix_lfc,
    sparseMatrix,
    package = "AcidTest",
    envir = environment()
)

df <- DFrame
gr <- GRanges
lfc <- matrix_lfc
mat <- matrix
rse <- RangedSummarizedExperiment
sce <- SingleCellExperiment
sparse <- sparseMatrix
txse <- SummarizedExperiment_transcripts

groceries <- c(NA, NA, "milk", "eggs", "eggs", "veggies")
mpgString <- "18.1, 18.7, 21, 21.4, 22.8"

organism <- "Homo sapiens"
release <- 97L

## nolint start
DataFrame <- S4Vectors::DataFrame
GRanges <- GenomicRanges::GRanges
IRanges <- IRanges::IRanges
SingleCellExperiment <- SingleCellExperiment::SingleCellExperiment
SummarizedExperiment <- SummarizedExperiment::SummarizedExperiment
`assay<-` <- SummarizedExperiment::`assay<-`
cause <- goalie::cause
hasInternet <- goalie::hasInternet
rowRanges <- SummarizedExperiment::rowRanges
skip_on_docker <- goalie::skip_on_docker
str_pad <- stringr::str_pad
tibble <- tibble::tibble
## nolint end

options(
    acid.save.dir = ".",
    acid.save.ext = "rds",
    acid.test = TRUE
)
