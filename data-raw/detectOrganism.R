## nolint start
suppressPackageStartupMessages({
    library(pipette)
})
## nolint end
object <- import("detectOrganism.csv")
object <- as(object, "DFrame")
saveRDS(object = object, file = "detectOrganism.rds")
export(object, "detectOrganism.csv")
