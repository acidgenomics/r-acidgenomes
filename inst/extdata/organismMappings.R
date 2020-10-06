library(basejump)
organismMappings <- as_tibble(import("organismMappings.csv"))
saveData(organismMappings, dir = ".")
