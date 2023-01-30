## https://www.gencodegenes.org/human/releases.html
## https://www.gencodegenes.org/mouse/releases.html
## nolint start
suppressPackageStartupMessages({
    library(pipette)
})
## nolint end
object <- import("gencode-to-ensembl.csv")
object <- as(object, "DataFrame")
saveRDS(object = object, file = "gencode-to-ensembl.rds")
