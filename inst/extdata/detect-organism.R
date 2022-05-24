## nolint start
suppressPackageStartupMessages({
    library(pipette)
})
## nolint end
object <- import("detect-organism.csv")
object <- as(object, "DataFrame")
saveRDS(object = object, file = "detect-organism.rds")
