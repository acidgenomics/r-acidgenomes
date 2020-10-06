library(basejump)
ucscMatrix <- as_tibble(import("ucscMatrix.csv"))
saveData(ucscMatrix, dir = ".")
