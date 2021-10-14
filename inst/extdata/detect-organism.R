library(basejump)
object <- import("detect-organism.csv")
object <- as(object, "DFrame")
saveRDS(object = object, file = "detect-organism.rds")
