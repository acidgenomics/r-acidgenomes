library(basejump)
object <- import("map-genome-build.csv")
object <- as(object, "DataFrame")
saveRDS(object = object, file = "map-genome-build.rds")
