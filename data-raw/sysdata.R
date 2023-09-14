library(usethis)
sysdataDetectOrganism <- readRDS("detectOrganism.rds")
sysdataMapNcbiTaxId <- readRDS("mapNcbiTaxId.rds")
use_data(
    sysdataDetectOrganism,
    sysdataMapNcbiTaxId,
    internal = TRUE,
    overwrite = TRUE
)
