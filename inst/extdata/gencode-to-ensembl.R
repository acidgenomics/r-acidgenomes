## https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/_README.TXT
## nolint start
suppressPackageStartupMessages({
    library(pipette)
})
## nolint end
object <- import("gencode-to-ensembl.csv")
object <- as(object, "DataFrame")
saveRDS(object = object, file = "gencode-to-ensembl.rds")
