library(testthat)
library(AcidGenomes)
test_check("AcidGenomes")
## Run the resource intensive unit tests less frequently (e.g. per month).
## Currently, I'm setting this in `~/.Rsecrets` for supported machines.
if (isTRUE(getOption("acid.test.extra"))) {
    ## Disable our flag for quick CI checks (see `helper_globals.R`).
    options(acid.test = FALSE)
    test_dir(file.path("tests", "testthat-extra"))
}
