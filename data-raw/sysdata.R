## nolint start
suppressPackageStartupMessages({
    library(usethis)
})
## nolint end
sysdataDetectOrganism <- readRDS("detectOrganism.rds")
sysdataMapNcbiTaxId <- readRDS("mapNcbiTaxId.rds")
use_data(
    sysdataDetectOrganism,
    sysdataMapNcbiTaxId,
    internal = TRUE,
    overwrite = TRUE
)
